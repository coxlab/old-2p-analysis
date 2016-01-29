%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter file for 2p data analysis
%%% Personal copy should be made per researcher
%%% Specify data folder to work from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameters=set_parameters(exp_name)
% Example:
% exp_name='2015-08-07_AH03/FOV01';
% parameters=set_parameters(exp_name)

%%% Get info about the environment
[system.host,system.user,system.platform]=get_ID();


%%% Set general parameters that define how the rest of pipeline will behave
ROI_params.ROI_definition_nr=1;

switch system.host
    case 'dixie'
        switch system.user
            case 'benvermaercke'
                root_folder='/nas/volume1/2photon';
                ROI_params.ROI_definition_nr=2; 
            case ''
            otherwise
                error('Username unknown')
        end
    case 'BKRUNCH'
        switch system.user
            case 'ben'
                root_folder='Q:\data\'; % quadraraid
                ROI_params.ROI_definition_nr=2;
            case ''
            otherwise
                error('Username unknown')
        end    
    case 'COXLAB-2P'
        switch system.user
            case 'labuser'
                root_folder='D:\Dropbox (coxlab)'; % seagate external drive
            otherwise
                error('Username unknown')
        end
    case 'BEN-PC'
        switch system.user
            case 'LBP'
                root_folder='C:\Users\LBP\Dropbox (coxlab)';
        end
        
        
    case 'Bens-MacBook-Pro.local'
        switch system.user
            case 'benvermaercke'
                root_folder='/Users/benvermaercke/Dropbox (coxlab)/';
                ROI_params.ROI_definition_nr=2;
            otherwise
                error('Username unknown')
        end
        
    %case ''
    otherwise        
        error('Hostname unknown...')
end

%%% Set folder structure
dirs.root_folder=root_folder;
dirs.raw_folder=fullfile(root_folder,'2p-data/raw',exp_name);
dirs.rec_folder=fullfile(root_folder,'2p-data/rec',exp_name);
dirs.reg_folder=fullfile(root_folder,'2p-data/reg',exp_name);
dirs.MWorks_folder=fullfile(root_folder,'MWorks',exp_name);
dirs.eyeTracker_folder=fullfile(root_folder,'eyetracker',exp_name);
dirs.analysis_folder=fullfile(root_folder,'analysis',exp_name);

if isdir(dirs.raw_folder)==0
    warning('Data folder not found on this machine...') 
else
    files=scandir(dirs.raw_folder,'.tif');        
    fprintf('Folder found and contains %d files! \n',length(files))
end

%%% get script folder
dirs.script_dir=fileparts(mfilename('fullpath'));
addpath(genpath(dirs.script_dir))

%%% Set popular look-up tables
colormaps.red=[linspace(0,1,256)' zeros(256,1) zeros(256,1)];
colormaps.green=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];

%%% Construct output
parameters.system=system;
parameters.dirs=dirs;
parameters.ROI_params=ROI_params;
parameters.colormaps=colormaps;



