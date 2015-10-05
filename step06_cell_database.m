clear all
clc

header_script
animal_ID='AH03';

offset_correction=[0 0];
switch animal_ID
    case 'AH02'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH02_20150803.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH02_resaved_im.png';
    case 'AH03'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH03_20150807.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH03_im.png';
        offset_correction=[-.4 .5];
    case 'AH05'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH05_20150901.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-09-01_AH05_im.png';
end

dataset_folder=fullfile(dataset_root,animal_ID);

dataset_files=scandir(dataset_folder,'mat');

for iFile=3%1:nFiles
    load_name=fullfile(dataset_folder,dataset_files(iFile).name);
    if exist(load_name,'file')
        load(load_name,'dataset')
        
        var=4;
        ROI_definitions=dataset.ROI_definitions;
        nROIs=length(ROI_definitions);
        for iROI=1:nROIs
            
            %%% Create a cell object based on the cell_data class, to hold
            %%% all data and perform all subsequent analyses.
            cell_data(iROI)=cell_processor(iROI,dataset);
            %%% Build condition matrix
            cell_data(iROI).build_condition_matrix(dataset.STIM)
            trace=dataset.RESP(:,iROI);
            cell_data(iROI).add_trace(trace)
            
            cell_data(iROI).do_RF_analysis()
            cell_data(iROI).do_RF_analysis_shuffled()
            %%
            cell_data(iROI).do_threshold(1.5)
            
            %cell_data(iROI).show_RF_map(cell_data(iROI).RF_map_TH)
            
            %%            
            cell_data(iROI).do_stimSelect_analysis()
            cell_data(iROI).get_sparseness()
            cell_data(iROI).calc_invariance()
           
            
            %             %% apply randomization procedure to account for noise
            %             nPerm=100;
            %             for iPerm=1:nPerm
            %                 random_trace=trace(randperm(length(trace)));
            %
            %                 trials_shuffled=struct;
            %                 for iTrial=1:nTrials
            %                     trial_data=condition_matrix(iTrial,:);
            %
            %                     %trace(trial_data(2):trial_data(3))
            %
            %                     trials_shuffled(iTrial).trial_nr=iTrial;
            %                     trials_shuffled(iTrial).condition_nr=trial_data(5);
            %                     trials_shuffled(iTrial).stim_nr=trial_data(6);
            %                     trials_shuffled(iTrial).pos_nr=trial_data(7);
            %                     trials_shuffled(iTrial).frames=trial_data(2)+1:trial_data(3)+3;
            %                     trials_shuffled(iTrial).nFrames=length(trials_shuffled(iTrial).frames);
            %                     trials_shuffled(iTrial).response=random_trace(trials_shuffled(iTrial).frames);
            %                     trials_shuffled(iTrial).response_avg=mean(trials_shuffled(iTrial).response);
            %                     trials_shuffled(iTrial).response_std=std(trials_shuffled(iTrial).response);
            %                 end
            %                 B=[cat(1,trials_shuffled.trial_nr) cat(1,trials_shuffled.condition_nr) cat(1,trials_shuffled.stim_nr) cat(1,trials.pos_nr) cat(1,trials_shuffled.response_avg) cat(1,trials_shuffled.response_std)];
            %                 results=[pivotTable(B,var,'mean',var) pivotTable(B,var,'mean',5) pivotTable(B,var,'std',5)];
            %                 RF_shuffled=reshape(results(:,2),4,8);
            %
            %                 mu_vector(iPerm)=mean(RF_shuffled(:));
            %                 sigma_vector(iPerm)=std(RF_shuffled(:));
            %             end
            %
            %             %%
            %             RF_norm=(RF-mean(mu_vector))/mean(sigma_vector);
            %             TH=1;
            %             RF_thresh=RF>TH;
            
            %             %% find best position
            %             [m,best_pos]=max(RF(:))
            %
            %             %% check stim selectivity there
            %             var=3;
            %             S=A(A(:,4)==best_pos,:);
            %             results=[pivotTable(S,var,'mean',var) pivotTable(S,var,'length',var) pivotTable(S,var,'mean',5) pivotTable(S,var,'std',5)];
            %
            %             %% find best stimulus
            %             [m,best_stim]=max(results(:,3))
            %             [sorted,order]=sort(results(:,3),'descend')
            %             %best_stim=[2]
            %             best_stim=order(1);
            %
            %             %%
            %             P=A(ismember(A(:,3),best_stim),:);
            %             var=4;
            %             results=[pivotTable(P,var,'mean',var) pivotTable(P,var,'length',var) pivotTable(P,var,'mean',5) pivotTable(P,var,'std',5)];
            %             results_padded=zeros(32,4);
            %             for iPos=1:32
            %                 sel=results(:,1)==iPos;
            %                 if any(sel)
            %                     results_padded(iPos,:)=results(sel,:);
            %                 else
            %                     results_padded(iPos,:)=[iPos 0 0 0];
            %                 end
            %             end
            %
            %             RF_best_stim=reshape(results_padded(:,3),4,8);
            %
            %             %% store data
            %
            %             cell_data.trace=trace;
            %             cell_data.condition_matrix=condition_matrix;
            %             cell_data.trials=trials;
            %             %%
            %             %figure
            %             %bar(A(:,5))
            %             %axis([0 11 -4 5])
            %
            %             imagesc(RF_best_stim)
            %             axis xy
            %
            %             R=sort(RF(:),'descend');
            %             %plot(R)
            %             set(gca,'cLim',[-.5 2])
            %
        end
        
        if 0
            %%
            subplot(121)
            INV=cat(1,cell_data.invariance_avg);
            hist(INV)
            
            %%
            subplot(122)
            RF_all=cat(3,cell_data.RF_map_TH);
            imagesc(mean(RF_all,3))
            colormap parula
            axis equal
            axis tight
            axis xy
        end
    end
end
