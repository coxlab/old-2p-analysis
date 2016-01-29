clear all
clc

%%% Set a folder containing raw tiffs you want to crunch
exp_name='2015-08-07_AH03/FOV01';

%%% Get parameters for running the rest of the pipeline
parameters=set_parameters(exp_name);

%%% Split data and flyback line, dump in respective folders (rec and MWorks)
%%% for res data: reconstruct

%%% Do motion correction using parameters set earlier

%%% Automatically detect ROIs, using selected strategy

%%% Manual check

%%% Extract timeseries data from registered images

%%% Combine with stimulus and behavioral data
