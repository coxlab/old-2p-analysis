clear all
clc

header_script

dataset_folder=fullfile(dataset_root,animal_ID);

files=scandir(data_folder,'tif');
nFiles=length(files);
if nFiles==0
    files=scandir(fullfile(data_folder,'data_analysis'),'mat');
    nFiles=length(files);
end

%%

file_name=fullfile(data_folder,files(1).name);
fprintf('Joining datasets in folder of file %s...\n',file_name)
save_name=fullfile(data_folder,'data_analysis',files(1).name);
save_name=strrep(save_name,'tif','mat');
load(save_name,'session_data')

%%% Make sure filenames are relative to data folder on this machine
session_data.rebase(data_root)
save_folder=session_data.folder_info.save_folder;

%% Join datasets

% check folder for joinable sessions: close in space (and time)
[clusters,cluster_vector,nClusters,files]=session_data.find_FOV_clusters(session_data.folder_info.save_folder);
%clusters([2 5])=NaN
%[(1:length(clusters))' clusters]
%die

% if close in space, but not in time, link files but call it a different
% timepoint, so you can look at evolution and/or stability

valid_datasets=zeros(nClusters,1);
for iClust=1:nClusters
    fprintf('Processing FOV %d\n',iClust)
    try
        cluster_nr=cluster_vector(iClust);
        %cluster_nr=clusters(iClust);
        session_vector=find(clusters==cluster_nr);
        nSessions=length(session_vector);
        
        %% Load all session into one compound object
        clear S
        for iSession=1:nSessions
            session_nr=session_vector(iSession);
            load_name=fullfile(save_folder,files(session_nr).name);            
            load(load_name,'session_data')
            session_data.rebase(data_root)
            session_data.ROI_definition_nr=ROI_definition_nr;
            S(iSession)=session_data;
        end
        disp('> Succeeded loading data sessions')
        
        
        %% Check whether same number of ROIs are defined in all sessions
        try
            S.get_ROI_counts(ROI_definition_nr)
            ROI_counts=S.get_ROI_counts(ROI_definition_nr);
        catch
            rethrow(lasterror)
            %error('error getting ROI counts...')
        end
        if all(ROI_counts==ROI_counts(1))
            nROIs=mean(ROI_counts);
        else
            ROI_counts
            error('Number of ROIs has to be identical for all sessions');
        end
        fprintf('> ROI numbers match: %d\n',nROIs)
        
        %% Verify distances are closely matched over sessions
        distances=S.get_ROI_distances(ROI_definition_nr);
        if any(mean(distances)>10)
            error('Large shift in ROI centers detected...')
        else
            disp('> ROI overlap checked completed successfully')
        end
        
        %% Verify that stim_matrix is healthy
        for iSession=1:nSessions
            if isfield(S(iSession).Experiment_info,'stim_matrix')
                if isempty(S(iSession).Experiment_info.stim_matrix)
                	error('Stim data is empty. Check step04...')
                end
            else
                error('No stim data found! Check step04...') 
            end
        end
        disp('> Stim data is healthy')
        
        %% Check if activity matrix matches ROI number
        for iSession=1:nSessions
            if size(S(iSession).Activity_traces.activity_matrix,2)~=nROIs
                error('Number of extracted traces does not match number of defined ROIs, run step04 again...')
            end
        end
        disp('> Number of extracted traces matches the number of ROIs')
        
        %% Then join them
        exp_name_str=exp_name;
        exp_name_str=strrep(exp_name_str,'/','_');
        
        dataset=S.join_data_sessions();
        
        %%% Add more features
        dataset.session_name=exp_name_str;
        dataset.animal_ID=animal_ID;        
        dataset.session_date=S(1).mov_info.mov_start_time; % use datevec to convert to numbers
        dataset.cluster_nr=iClust;
        dataset.session_vector=session_vector;
        dataset.nSessions=length(session_vector);
        dataset.nFrames=size(dataset.STIM,1);
        dataset.frame_rate=S(1).mov_info.frame_rate;
        %valid_datasets(iClust)=1;
        
        if save_it==1            
            save_name=fullfile(dataset_folder,sprintf([exp_name_str '_FOV%02d.mat'],iClust));
            savec(save_name)
            save(save_name,'dataset')
            disp('> Saving the dataset')
        end
        
        
        
    catch
        %%
        rethrow(lasterror)
        A=lasterror;
        disp('Unable to create dataset from these session files...')
        disp(A.message)
        session_vector
    end
end



%%% Save database file and put all subsequent code into imaging_datasets
%%% class.



%
%
% %%% Do something with dataset
% valid_dataset_vector=find(valid_datasets==1);
%
% if 0
%     %% Plot FOV, works on 1 or more datasets
%     dataset(valid_dataset_vector).plot_FOV()
% end
% if 0
%     %% look at temporal profile
%     dataset_nr=2;
%     dataset(dataset_nr).plot_traces()
% end
%
% if 0
%     %% Look at RF of individual ROI
%     dataset_nr=2;
%     ROI_nr=4;
%     %dataset(dataset_nr).RF_analysis(ROI_nr)
%     dataset(dataset_nr).plot_RF_map(ROI_nr,0)
% end
%
% %% Plot RF maps for all ROIs
% N=length(valid_dataset_vector);
% IND_AVG=2;
% deconvolve=0;
% for iDS=1:N
%     dataset_nr=valid_dataset_vector(iDS);
%
%     switch IND_AVG
%         case 1
%             %%% Show individual RF maps
%             dataset(dataset_nr).plot_RF_map([],deconvolve)
%         case 2
%             %%% Show average map
%             results=dataset(dataset_nr).RF_analysis([],deconvolve);
%             AVG_MAP=mean(cat(3,results.condAverage.RF_map),3);
%             AVG_MAP=imresize(AVG_MAP,10);
%
%             figure(dataset_nr)
%             clf
%             imagesc(AVG_MAP)
%             axis xy
%             title(sprintf('Dataset #%d',dataset_nr))
%     end
% end
%
%
% %% Think about how to show RF of each neuron in 1 representation, or 2.
% % Give azimuth 1 color on a map and then elevation on another.
% % 1 color point per neuron
% % Apply some kind of smoothing to get coherent maps (interp?)
%
% % need chance level to use a fixed TH
% % then calc weighted average of above TH point to get 'center' of RF
% % use x-y position of those centers to create azimuth and elevation maps
% figure()
% dataset_nr=3;
% TH=0.5;
%
% data_matrix=zeros(0,5);
% count=1;
% for iROI=1:dataset(dataset_nr).nROIs
%     ROI=dataset(dataset_nr).ROI_definitions(iROI);
%     results=dataset(dataset_nr).RF_analysis(iROI,0);
%     map=results.condAverage.RF_map;
%     map_TH=map>TH;
%     %map_TH=imopen(map_TH,strel('disk',1));
%
%     [x,y]=find(map_TH);
%     N=length(x);
%     if N>1
%         values=zeros(N,1);
%         for iRow=1:N
%             values(iRow)=map(x(iRow),y(iRow));
%         end
%         M=[x y values/sum(values)];
%         center_Y=sum(M(:,1).*M(:,3));
%         center_X=sum(M(:,2).*M(:,3));
%
%         data_matrix(count,:)=[iROI ROI.center_coords center_X center_Y]
%         count=count+1;
%
%         imagesc(map_TH)
%         hold on
%         plot(center_X,center_Y,'rs')
%         hold off
%         axis xy
%     else
%         % untuned
%         clf
%     end
% end
%
%
% %%% AZIMUTH= col 4
% nConditions=8;
% nPoints=size(data_matrix,1);
% all_coords=cat(1,dataset(dataset_nr).ROI_definitions.center_coords);
%
% color_range=jet(80);
%
% subplot(2,3,[1 3])
% imshow(calc_gamma(dataset(dataset_nr).MIP_std.data,.3),[])
% hold on
% plot(all_coords(:,1),all_coords(:,2),'r.')
%
% for iPoint=1:nPoints
%     loc=round(data_matrix(iPoint,4)*10);
%     plot(data_matrix(iPoint,2),data_matrix(iPoint,3),'o','color',color_range(loc,:))
% end
% colormap(green)
% colorbar
% %set(gca,'CLim',[1 nConditions])
% hold off
% axis ij
%
% subplot(2,3,4)
% hist(data_matrix(:,4))
% axis([1 8 0 4])
% subplot(2,3,5)
% hist(data_matrix(:,5))
% axis([1 4 0 4])
% subplot(2,3,6)
% plot(data_matrix(:,4),data_matrix(:,5),'rx')
% axis([1 8 1 4])
% axis equal
% axis square
%
%
