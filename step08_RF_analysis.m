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
condition_vector=stim_matrix(stim_matrix(:,5)>-1,col_nr);
unique_conditions=unique(condition_vector);
nConditions=length(unique_conditions);

%%% For each ROI, determine response per condition
stim_length_frames=9;
cond_matrix_mean=zeros(nConditions,nROI);
cond_matrix_std=cond_matrix_mean;
perm_matrix_mean=cond_matrix_mean;
perm_matrix_std=cond_matrix_mean;
for iROI=1:nROI
    for iCond=1:nConditions
        condition_nr=unique_conditions(iCond);
        
        %%% Real data
        sel=stim_matrix(:,col_nr)==condition_nr;
        
        %%% find different presentation of condition
        frames_selected=find(sel);
        start_frames=[0;diff(frames_selected)>1];
        repeat_starts=frames_selected(start_frames==1);
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
        
        %%% Perm data: do multiple times, needs reformatting of the loop
        stim_matrix_shuffled=Shuffle(stim_matrix(:,col_nr));
        sel_rand=stim_matrix_shuffled==condition_nr;
        
        %%% find different presentation of condition
        frames_selected=find(sel_rand);
        start_frames=[0;diff(frames_selected)>1];
        repeat_starts=frames_selected(start_frames==1);
        repeat_starts(repeat_starts+stim_length_frames>nFrames)=[];
        
        nRepeats=length(repeat_starts);
        repeat_vector=zeros(nRepeats,1);
        for iRepeat=1:nRepeats
            repeat_start=repeat_starts(iRepeat);
            repeat_vector(iRepeat)=mean(resp_matrix(repeat_start:repeat_start+7,iROI));
            %stim_matrix(repeat_start-2:repeat_start+6,:); % sanity check
        end
        perm_matrix_mean(iCond,iROI)=mean(repeat_vector);
        perm_matrix_std(iCond,iROI)=std(repeat_vector);
    end
end

%%
clf
nCols=ceil(sqrt(nROI));
nRows=ceil(nROI/nCols);

MIP=dataset.MIP.data;

for iROI=1:nROI
    subplot(nRows,nCols,iROI)
    switch col_nr
        case 5
            RF=flipud(reshape(cond_matrix_mean(:,iROI),3,4));
        case 8            
            RF=flipud(reshape(cond_matrix_mean(:,iROI),4,8));
    end
    mu=mean(perm_matrix_mean(:,iROI));
    sigma=std(perm_matrix_mean(:,iROI));
    RF_zscore=(RF-mu)/sigma;
    
    %%% Place RF map in FOV    
    center=round(dataset.ROI_definitions(iROI).center_coords);
    RF_disp=imresize(RF_zscore*700,3,'bicubic');
    MIP(center(2)-size(RF_disp,1)/2+1:center(2)+size(RF_disp,1)/2,center(1)-size(RF_disp,2)/2+1:center(1)+size(RF_disp,2)/2)=RF_disp;
    
    switch 2
        case 1
            imagesc(RF_zscore)
            set(gca,'Clim',[-1 1]*25)
        case 2
            imagesc(RF_zscore>10)
            set(gca,'Clim',[0 1])
    end
    title(dataset.ROI_definitions(iROI).center_coords)
    axis equal
    axis tight
    set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
end
clf
imshow(MIP,[])
colormap(green)


