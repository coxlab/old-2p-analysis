clear all
clc

header_script
load(fullfile(data_folder,'data_analysis','2015-08-06_AH02_001.mat'))

%%
[CC_max,offset]=im_align(session_data.MIP_avg.data,session_data.MIP_cc_local.data);

%%
[CC_max,offset]=im_align(session_data.MIP_avg.data,session_data.MIP_cc_local.data,1);