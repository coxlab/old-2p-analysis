%%% When on mac, working off of dropbox
%% BV20150420: when on linux server, use different
if isunix==1
    if ismac
        [~, user_name] = system('whoami');user_name=user_name(1:end-1);
        root_folder=fullfile('/Users',user_name,'/Dropbox (coxlab)');
    else
        root_folder='/share/scratch2/benv/2p_data';
    end
else % 2p
    [~, user_name] = system('whoami');user_name=user_name(1:end-1);
    root_folder=fullfile('/Users',user_name,'/Dropbox (coxlab)');
end

% define data folder
data_root=fullfile(root_folder,'2p-data');
mworks_folder=fullfile(root_folder,'Analysis/Scripting/Matlab');

% Add code for MWK analysis scripts
addpath(genpath(mworks_folder))
path_dir=fileparts(mfilename('fullpath'));
addpath(genpath(path_dir))

green=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];

