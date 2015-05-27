clear all
clc

%%% for the most responsive position, do selectivity analysis over each
%%% cell or each population. Try also to see how generalization works over
%%% all other positions

%%% BV20150518: unbiased method: every 
% find number of presentations per unique combination of shape and position
% i.e. 1/384 and concat all the responses per FOV 

header_script
save_it=0;
dataset_selector=1;

%%% Load requested merged dataset
loadName=fullfile(data_folder,'data_analysis',sprintf('dataset_%03d.mat',dataset_selector));
load(loadName,'dataset')

nFrames=dataset.nFrames
nROI=dataset.nROI;
stim_matrix=dataset.stim_matrix;
resp_matrix=dataset.resp_matrix;
resp_matrix_NND=dataset.resp_matrix_NND;


%%
stim_matrix_crop=stim_matrix;
stim_matrix_crop(stim_matrix_crop(:,4)==-1,:)=[];
M=[pivotTable(stim_matrix_crop,8,'mean',8) pivotTable(stim_matrix_crop,8,'mean',5) pivotTable(stim_matrix_crop,8,'length',8)];


M
size(M)

%% Collect average responses for all ROIs



%%




