clear all
clc

header_script

files=scandir(data_folder,'tif');
nFiles=length(files);
%%

file_name=fullfile(data_folder,files(1).name);
fprintf('Joining datasets in folder of file %s...\n',file_name)
save_name=fullfile(data_folder,'data_analysis',files(1).name);
save_name=strrep(save_name,'tif','mat');
load(save_name,'session_data')

%%% Make sure filenames are relative to data folder on this machine
session_data.rebase(data_root)


%% Join datasets

% check folder for joinable sessions: close in space (and time)
[clusters,cluster_vector,nClusters,files]=session_data.find_FOV_clusters(session_data.folder_info.save_folder);
clusters(22)=NaN

% if close in space, but not in time, link files but call it a different
% timepoint, so you can look at evolution and/or stability

valid_datasets=zeros(nClusters,1);
for iClust=1:nClusters
    try
        cluster_nr=cluster_vector(iClust);
        session_vector=find(clusters==cluster_nr);
        nSessions=length(session_vector);
        
        %% Load all session into one compound object
        clear S
        for iSession=1:nSessions
            session_nr=session_vector(iSession);
            load_name=fullfile(session_data.folder_info.save_folder,files(session_nr).name);
            load(load_name,'session_data')
            S(iSession)=session_data;
        end
        
        %% Check whether same number of ROIs are defined in all sessions
        ROI_counts=S.get_ROI_counts(ROI_definition_nr);
        if all(ROI_counts==ROI_counts(1))
            nROIs=mean(ROI_counts);
        else
            error('Number of ROIs has to be identical for all sessions');
        end
        
        
        %% Verify distances are closely matched over sessions
        distances=S.get_ROI_distances(ROI_definition_nr);
        if any(mean(distances)>10)
            error('Large shift in ROI centers detected...')
        else
            % all good: take the ROI_definitions of one of the sessions
            %ROI_definitions=S(1).ROI_definitions(ROI_definition_nr).ROI;
            
            % Then join them
            % leave out last trial
            dataset(iClust)=S.join_data_sessions();
            dataset(iClust).cluster_nr=iClust;
            dataset(iClust).session_vector=session_vector;
            dataset(iClust).nSessions=length(session_vector);
            dataset(iClust).nFrames=size(dataset(iClust).STIM,1);
            dataset(iClust).frame_rate=S(1).mov_info.frame_rate;
            valid_datasets(iClust)=1;
        end
        
    catch
        disp('Unable to create dataset from these session files...')
        session_vector
    end
end



%%% Do something with dataset
valid_dataset_vector=find(valid_datasets==1);
N=length(valid_dataset_vector);


if 0
    %% Plot FOV, works on 1 or more datasets
    dataset(valid_dataset_vector).plot_FOV()
end
if 0
    %% look at temporal profile
    dataset_nr=2;
    dataset(dataset_nr).plot_traces()
end

if 0
    %% Look at RF of individual ROI
    dataset_nr=1;
    ROI_nr=36;
    %dataset(dataset_nr).RF_analysis(ROI_nr)
    dataset(dataset_nr).plot_RF_map(ROI_nr,0)
end

%% Plot RF maps for all ROIs
IND_AVG=2;
deconvolve=0;
for iDS=1:N
    dataset_nr=valid_dataset_vector(iDS);    
    
    switch IND_AVG
        case 1
            %%% Show individual RF maps
            dataset(dataset_nr).plot_RF_map([],deconvolve)
        case 2
            %%% Show average map
            results=dataset(dataset_nr).RF_analysis([],deconvolve);
            AVG_MAP=mean(cat(3,results.condAverage.RF_map),3);
            AVG_MAP=imresize(AVG_MAP,10);
            
            figure(dataset_nr)
            clf
            imagesc(AVG_MAP)
            axis xy
            title(sprintf('Dataset #%d',dataset_nr))
    end    
end

