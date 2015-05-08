clear all
clc

header_script

dataset_selector=2;

%%% Load requested merged dataset
loadName=fullfile(data_folder,'data_analysis',sprintf('dataset_%03d.mat',dataset_selector));
load(loadName,'dataset')

nFrames=dataset.nFrames;
nROI=dataset.nROI;
stim_matrix=dataset.stim_matrix;
resp_matrix=dataset.resp_matrix;
resp_matrix_NND=dataset.resp_matrix_NND;

%%% Get unique condition numbers
col_nr=8;
condition_vector=stim_matrix(stim_matrix(:,col_nr)>0,col_nr);
unique_conditions=unique(condition_vector);
nConditions=length(unique_conditions);

%%% For each ROI, determine response per condition
stim_length_frames=9;
cond_matrix_mean=zeros(nConditions,nROI);
cond_matrix_std=cond_matrix_mean;
for iROI=1:nROI
    for iCond=1:nConditions
        condition_nr=unique_conditions(iCond);
        sel=stim_matrix(:,col_nr)==condition_nr;
        
        %%% find different presentation of condition
        frames_selected=find(sel);
        start_frames=find([0;diff(frames_selected)>1]);
        repeat_starts=frames_selected(start_frames);
        repeat_starts(repeat_starts+stim_length_frames>nFrames)=[];
        
        nRepeats=length(repeat_starts);
        repeat_vector=zeros(nRepeats,1);
        for iRepeat=1:nRepeats
            repeat_start=repeat_starts(iRepeat);
            repeat_vector(iRepeat)=mean(resp_matrix(repeat_start+2:repeat_start+7,iROI));
            %stim_matrix(repeat_start-2:repeat_start+6,:); % sanity check
        end
        cond_matrix_mean(iCond,iROI)=mean(repeat_vector);
        cond_matrix_std(iCond,iROI)=std(repeat_vector);
    end
end

%%
nCols=ceil(sqrt(nROI));
nRows=ceil(nROI/nCols);
for iROI=1:nROI
    subplot(nRows,nCols,iROI)
    RF=flipud(reshape(cond_matrix_mean(:,iROI),4,8));
    imagesc(RF)
    set(gca,'Clim',[-1 1]*5)
    axis equal
    axis tight
end



