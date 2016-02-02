clear all
clc

%%%
% A 2-photon experiment usually exists of 1 or more raw data files (e.g.
% tif), (visual) stimulation records and multiple datastreams that are
% collected simultaneously. This package aims to provide a unifying
% approach to analyse this type of data.
% First, we will process, motion-correct movies and then select ROIs using
% a variety of methods. The extracted time traces of these ROIs will then
% be combined with stimulation and other datastreams.

%%% Set a folder containing raw tiffs you want to crunch
%exp_name='2015-08-07_AH03/FOV01';

exp_name='160201_KS159_2P_KS/run03_ori8_V1_random';
lab_name='boninlab';

%%% Get parameters for running the rest of the pipeline
%parameters=set_parameters(exp_name); % function
dataset=data_object(exp_name);   % object
%dataset.toggle_save();
%%

switch lab_name
    case 'coxlab'
        %%% For coxlab galvo data:
        %%% Split data and flyback line, 
        %split_data(parameters)
        dataset.preprocess_data()        
        
        %%% Dump in respective folders (rec and MWorks)
        %parameters.dirs.rec_folder
        %parameters.dirs.MWorks_folder        
        
    case 'boninlab'
        %%% For bonin lab res data: reconstruct
        dataset.reconstruct_data()
end
dataset.save_data()

%% Do motion correction using parameters set earlier
dataset.register_data()

%%% Automatically detect ROIs, using selected strategy

%%% Manual check

%%% Extract timeseries data from registered images

%%% Combine with stimulus and behavioral data
