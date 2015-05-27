%%% Export the result of this step so we can build on it in the next step
clear all
clc

header_script

%%
save_it=1;
rehash=0;

col_nr=8;
map_size=[4 8];

nPerm=100;
shuffle_method=2;

blow_fields=4;
blow_field_avg=50;
z_scale=[-1 1]*4;
z_scale_AVG=[-1 1]*.75;


for dataset_selector=1:3
    for data_type=1                        
        %%% Load requested merged dataset
        loadName=fullfile(data_folder,'data_analysis',sprintf('dataset_%03d.mat',dataset_selector));
        load(loadName,'dataset')
        
        if isfield(dataset,'MAPs')&&rehash==0
            disp('reloading data')
        else
            nFrames=dataset.nFrames;
            nROI=dataset.nROI;
            stim_matrix=dataset.stim_matrix;
            resp_matrix=dataset.resp_matrix;
            resp_matrix_NND=dataset.resp_matrix_NND;
            
            
            %%% Rewrite analysis procedure
            % Collect the input data, stim matrix and responses
            % Feed these to your analysis routine, dividing up per condition and
            % getting stats
            % Then shuffle the input and reanalyse a number of times
            
            %%% Define relevant condition
            groups=stim_matrix(:,col_nr); % position
            groups(groups==0)=-1; % do not do this for stim type! There 0 is stim01
            
            %%% Get raw maps
            MAPs=struct;
            for iROI=1:nROI
                % Get response trace
                switch data_type
                    case 1
                        Y=resp_matrix(:,iROI);
                    case 2
                        Y=resp_matrix_NND(:,iROI);
                end
                
                % Do calculations and return result
                result=analyse_RF(Y,groups);
                
                MAP=flipud(reshape(cat(1,result.mean_resp),map_size(1),map_size(2)));
                MAPs(iROI).MAP=MAP;
            end
            
            %%% Randomize and recalculate
            t0=clock;
            for iROI=1:nROI % 46 is a nice example
                MU_vector=zeros(nPerm,1);
                SIGMA_vector=MU_vector;
                for iPerm=1:nPerm
                    %% Shuffle stim matrix
                    switch shuffle_method
                        case 1 % naive shuffle
                            groups_shuffled=groups(randperm(nFrames),:);
                        case 2 % shuffle stim identity
                            groups_shuffled=smart_shuffle(groups);
                    end
                    
                    % Get response trace
                    switch data_type
                        case 1
                            Y=resp_matrix(:,iROI);
                        case 2
                            Y=resp_matrix_NND(:,iROI);
                    end
                    
                    % Do calculations and return result
                    result=analyse_RF(Y,groups_shuffled);
                    MU_vector(iPerm)=mean(cat(1,result.mean_resp));
                    SIGMA_vector(iPerm)=mean(cat(1,result.std_resp));
                    
                    if plot_it==1
                        %%% check intermediate results
                        MAP=flipud(reshape(cat(1,result.mean_resp),map_size(1),map_size(2)));
                        R=MAPs(:,:,iROI);
                        disp([std(R(:)) std(MAP(:))])
                        max_val=max([std(R(:)) std(MAP(:))])*3;
                        
                        subplot(311)
                        plot(Y)
                        subplot(312)
                        imshow(R,[-max_val max_val])
                        title(max_val)
                        subplot(313)
                        imshow(MAP,[-max_val max_val])
                    end
                end
                %std(MU_vector)
                MAPs(iROI).nPerm=nPerm;
                MAPs(iROI).MU=mean(MU_vector);
                MAPs(iROI).MU_std=std(MU_vector);
                MAPs(iROI).SIGMA=mean(SIGMA_vector);
                MAPs(iROI).SIGMA_std=std(SIGMA_vector);
                
                %%% Show progress
                progress(iROI,nROI,t0)
            end
            
            %%% Calc z-scored maps : weigh each RF map by how much it differs from noise
            for iROI=1:nROI
                MAP=MAPs(iROI);
                MAP_zscored=(MAP.MAP-MAP.MU)/MAP.SIGMA;
                MAPs(iROI).MAP_zscored=MAP_zscored;
            end
            
            switch data_type
                case 1
                    dataset.MAPs=MAPs;
                case 2
                    dataset.MAPs_NND=MAPs;
            end
            
            if save_it==1
                %% Save data to dataset file
                save(loadName,'dataset')
                disp(['Updated ' loadName])
            end
        end
        
        switch data_type
            case 1
                MAPs=dataset.MAPs;
            case 2
                MAPs=dataset.MAPs_NND;
        end
        
        
        %% show mip
        MIP_RGB=cat(3,dataset.MIP.data*0,dataset.MIP.data,dataset.MIP.data*0);
        MIP_RGB=calc_gamma(MIP_RGB,.5);
        MIP_RGB=MIP_RGB/max(MIP_RGB(:));
        
        imshow(MIP_RGB,[])
        axis([0 512 0 399])
        
        if save_it
            %%
            parts=strsplit(data_folder,filesep);
            saveName=fullfile(data_folder,'data_analysis','RF_maps',sprintf([parts{end} '_%03d_MIP.tif'],dataset_selector));
            savec(saveName)
            imwrite(MIP_RGB,saveName)
        end
        
        
        %% Show mean
        mean_RF=mean(cat(3,MAPs.MAP_zscored),3);
        mean_RF_blowup=imresize(mean_RF,blow_field_avg);
        %dataset.mean_RF=mean_RF;
        %dataset.mean_RF_blowup=mean_RF_blowup;
        
        %MEAN_RGB=applyColorMap2(mean_RF_blowup,'hot',0,z_scale_AVG(2));
        MEAN_RGB=applyColorMap2(mean_RF_blowup,'hot',0,range(mean_RF_blowup(:))*1.2);
        
        imshow(MEAN_RGB,[])
        
        colormap hot
        axis equal
        box off
        axis tight
        set(gca,'Clim',[0 255])
        
        if save_it
            %%
            parts=strsplit(data_folder,filesep);
            saveName=fullfile(data_folder,'data_analysis','RF_maps',sprintf([parts{end} '_%03d_MEAN.tif'],dataset_selector));
            imwrite(uint8(MEAN_RGB*256),saveName)
        end
        
        %% Put RFs in a spatial arrangement
        %MIP=NaN(size(dataset.MIP.data));
        MIP=zeros(size(dataset.MIP.data))-z_scale(2);
        for iROI=1:dataset.nROI
            center=round(dataset.ROI_definitions(iROI).center_coords);
            RF=MAPs(iROI).MAP_zscored;
            RF_disp=imresize(RF,blow_fields,'bicubic');
            
            MIP(center(2)-size(RF_disp,1)/2+1:center(2)+size(RF_disp,1)/2,center(1)-size(RF_disp,2)/2+1:center(1)+size(RF_disp,2)/2)=RF_disp;
        end
        
        MIP=MIP/max(MIP(:));
        MIP_RGB=applyColorMap2(MIP,'',0,1);
        imshow(MIP_RGB,z_scale)
        %colormap hot
        
        if save_it
            %%
            parts=strsplit(data_folder,filesep);
            saveName=fullfile(data_folder,'data_analysis','RF_maps',sprintf([parts{end} '_%03d_MAP.tif'],dataset_selector));
            imwrite(MIP_RGB,saveName)
        end
        
        
    end
end
% %%
%
% die
%
% %% Get unique condition numbers
% condition_vector=stim_matrix(stim_matrix(:,5)>-1,col_nr);
% unique_conditions=unique(condition_vector);
% nConditions=length(unique_conditions);
%
% %%% For each ROI, determine response per condition
% stim_length_frames=24;
% cond_matrix_mean=zeros(nConditions,nROI);
% cond_matrix_std=cond_matrix_mean;
% perm_matrix_mean=cond_matrix_mean;
% perm_matrix_std=cond_matrix_mean;
% for iROI=1:nROI
%     for iCond=1:nConditions
%         condition_nr=unique_conditions(iCond);
%
%         %%% Real data
%         sel=stim_matrix(:,col_nr)==condition_nr;
%
%         %%% find different presentation of condition
%         frames_selected=find(sel);
%         start_frames=[1;diff(frames_selected)>1];
%         repeat_starts=frames_selected(start_frames==1);
%         repeat_starts(repeat_starts+stim_length_frames>nFrames)=[];
%
%         nRepeats=length(repeat_starts);
%         if nRepeats==0
%             die
%         end
%         repeat_vector=zeros(nRepeats,1);
%         for iRepeat=1:nRepeats
%             repeat_start=repeat_starts(iRepeat);
%             repeat_vector(iRepeat)=mean(resp_matrix(repeat_start+2:repeat_start+7,iROI));
%             %stim_matrix(repeat_start-2:repeat_start+6,:); % sanity check
%         end
%         cond_matrix_mean(iCond,iROI)=mean(repeat_vector);
%         cond_matrix_std(iCond,iROI)=std(repeat_vector);
%
%         %%% Perm data: do multiple times, needs reformatting of the loop
%         stim_matrix_shuffled=Shuffle(stim_matrix(:,col_nr));
%         sel_rand=stim_matrix_shuffled==condition_nr;
%
%         %%% find different presentation of condition
%         frames_selected=find(sel_rand);
%         start_frames=[0;diff(frames_selected)>1];
%         repeat_starts=frames_selected(start_frames==1);
%         repeat_starts(repeat_starts+stim_length_frames>nFrames)=[];
%
%         nRepeats=length(repeat_starts);
%         repeat_vector=zeros(nRepeats,1);
%         for iRepeat=1:nRepeats
%             repeat_start=repeat_starts(iRepeat);
%             repeat_vector(iRepeat)=mean(resp_matrix(repeat_start:repeat_start+7,iROI));
%             %stim_matrix(repeat_start-2:repeat_start+6,:); % sanity check
%         end
%         perm_matrix_mean(iCond,iROI)=mean(repeat_vector);
%         perm_matrix_std(iCond,iROI)=std(repeat_vector);
%     end
% end
%
% %%%
% clf
% nCols=ceil(sqrt(nROI));
% nRows=ceil(nROI/nCols);
%
% MIP=dataset.MIP.data*NaN;
%
% for iROI=1:nROI
%
%     switch col_nr
%         case 5
%             RF=flipud(reshape(cond_matrix_mean(:,iROI),3,4));
%         case 8
%             RF=flipud(reshape(cond_matrix_mean(:,iROI),4,8));
%     end
%     mu=mean(perm_matrix_mean(:,iROI));
%     sigma=std(perm_matrix_mean(:,iROI));
%     RF_zscore=(RF-mu)/sigma;
%
%     RF_zscore_all(:,:,iROI)=RF_zscore;
%
%     RF_zscore_disp=double(RF_zscore)/5*double(intmax('uint16'))/2+double(intmax('uint16'))/2;
%     %%% Place RF map in FOV
%     center=round(dataset.ROI_definitions(iROI).center_coords);
%     RF_disp=imresize(RF_zscore_disp,3,'bicubic');
%     MIP(center(2)-size(RF_disp,1)/2+1:center(2)+size(RF_disp,1)/2,center(1)-size(RF_disp,2)/2+1:center(1)+size(RF_disp,2)/2)=RF_disp;
%
%     if plot_it==1
%         subplot(nRows,nCols,iROI)
%         switch 2
%             case 1
%                 imagesc(RF_zscore)
%                 set(gca,'Clim',[-1 1]*25)
%             case 2
%                 imagesc(RF_zscore>10)
%                 set(gca,'Clim',[0 1])
%         end
%         title(dataset.ROI_definitions(iROI).center_coords)
%         axis equal
%         axis tight
%         set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
%     end
% end
%
% %%
% figure(dataset_selector)
% subplot(211)
% MIP_RGB=cat(3,dataset.MIP.data*0,dataset.MIP.data,dataset.MIP.data*0);
% MIP_RGB=calc_gamma(MIP_RGB,.5);
% MIP_RGB=MIP_RGB/max(MIP_RGB(:));
%
% imshow(MIP_RGB,[])
% %colormap(green)
% subplot(212)
% imshow(MIP,[])
% colormap(hot)
%
% if save_it==1
%     %%
%     parts=strsplit(data_folder,filesep);
%     saveName=fullfile(data_folder,'data_analysis','RF_maps',sprintf([parts{end} '_%03d.eps'],dataset_selector));
%     savec(saveName)
%     print(gcf,saveName,'-depsc')
% end
%
% if 0
%     %% show mean map
%     imagesc(mean(RF_zscore_all,3))
%     set(gca,'Clim',[-1 1]*5)
% end
% %% plot x and y centers values at cell location and interpolate between them, this should give you a smooth RF map
% % use color to code position
% % make one plot for each coord: azimuth and elevation
%
%
%
