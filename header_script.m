%%% When on mac, working off of dropbox
%% BV20150420: when on linux server, use different
if isunix==1
    if ismac
        [~, user_name] = system('whoami');user_name=user_name(1:end-1);
        root_folder=fullfile('/Users',user_name,'/Dropbox (coxlab)');
    else % on linux, define path to current location of dropbox folders
        [~, user_name] = system('whoami');user_name=user_name(1:end-1);
        root_folder=fullfile('/home/',user_name,'/Dropbox (coxlab)/'); % temp location
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

% either set data_folder in code or use uigetdir
switch 0
    case 0 % hardcoded
        % Ben data
        %exp_name='2015-03-03_AF03-light_awake';
        %exp_name='2015-03-04_AF11';
        %exp_name='2015-03-05_AF03';
        
        exp_name='2015-03-05_AF11';
        %exp_name='2015-04-07_AF11';
        
        %exp_name='2015-04-10_AF11_exp';
        %exp_name='2015-04-15_AF11';
        %exp_name='2015-04-16_AF11';
        
        
        % Julia data
        %exp_name='20150502_jrat3/Session02';
        %exp_name='20150502_jrat3/Session03-05';
        %exp_name='20150502_jrat3/Session06-10';
        %exp_name='20150502_jrat3/Session11-16';
        %exp_name='20150502_jrat3/Session17-22';
        %exp_name='20150502_jrat3/Session23-25';
        
        %data_folder=['/Users/' user_name '/Dropbox (coxlab)/2p-data/' exp_name];
        data_folder=fullfile(data_root,exp_name);
        if isdir(data_folder)==0
            data_folder
            error('directory not found')
        end
    case 1 % using uigetdir
        curr_dir=pwd;
        cd(data_root)
        data_folder=uigetdir(data_root);
end
cd(root_folder)

% suppress output so we can run in -nojvm mode
plot_it=0;
save_it=1;
use_custom_ROI_defs=1; % Relevant for step 3 (GUI)


% display current working folder
fprintf('Current working folder is : %s\n',data_folder)
