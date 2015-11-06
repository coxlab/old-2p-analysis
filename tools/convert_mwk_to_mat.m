%%% this script will convert .mwk data files to .mat files, everything or parts of it
%%% only works on mac systems so far.
%%% perhaps, we could build a python version that is portable across
%%% systems if python can export to .mat files

%%% in principle, we could crawl the entire data directory and convert
%%% every .mwk file we encounter. perhaps this will be the default behavior
%%% if just folder en not file is provided as input.

%%% We are in the process of trying to compile the mwfeval function for
%%% windows and linux, no luck so far... code seems to be outdated

clear all
clc

if ismac
    
else
    disp('Currently, this script only run on mac systems given that the mex file needed to read the .mwk file is only compiled for mac')
end