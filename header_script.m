%%% When on mac, working off of dropbox
%% BV20150420: when on linux server, use different
if isunix==1
    if ismac
        [~, user_name] = system('whoami');user_name=user_name(1:end-1);
        root_folder=fullfile('/Users',user_name,'/Dropbox (coxlab)');
    else % on linux, define path to current location of dropbox folders
        [~, user_name] = system('whoami');user_name=user_name(1:end-1);
        %root_folder=fullfile('/home/',user_name,'/Dropbox (coxlab)/'); % temp location
        root_folder='/nas/volume1/2photon/2p-data'; % temp location
    end
else % 2p windows
    user_name=getenv('username');
    switch getenv('computername')
        case 'BKRUNCH'
             root_folder='Q:\data\'; % quadraraid
        case 'BEN-PC'
            root_folder='C:\Users\LBP\Dropbox (coxlab)';
        case 'COXLAB-2P'
             root_folder='D:\Dropbox (coxlab)'; % seagate external drive
        otherwise
            [~, user_name] = system('whoami');user_name=user_name(1:end-1);
            switch user_name
                case 'ben-pc\lbp'
                    root_folder='C:\Users\LBP\Dropbox (coxlab)';
                otherwise
                    %root_folder=fullfile('/Users',user_name,'/Dropbox (coxlab)');
                    root_folder='D:\Dropbox (coxlab)'; % seagate external drive
            end
    end
end

% define data folder
switch getenv('computername')
    case 'BKRUNCH'
        data_root=fullfile(root_folder,'2photon\reg');
        mworks_folder=fullfile(root_folder,'Analysis/Scripting/Matlab');
        dataset_root=fullfile(root_folder,'2p-datasets');
    otherwise
        data_root=fullfile(root_folder,'2p-data');
        mworks_folder=fullfile(root_folder,'Analysis/Scripting/Matlab');
        dataset_root=fullfile(root_folder,'2p-datasets');
end

% Add code for MWK analysis scripts
addpath(genpath(mworks_folder))
path_dir=fileparts(mfilename('fullpath'));
addpath(genpath(path_dir))


tools_dir='/Users/benvermaercke/Bonin/2p_repo/tools';
addpath(genpath(tools_dir))
tools_dir='/Users/benvermaercke/Dropbox (coxlab)/tools';
addpath(genpath(tools_dir))

red=[linspace(0,1,256)' zeros(256,1) zeros(256,1)];
green=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];


% either set data_folder in code or use uigetdir
switch 0
    case 0 % hardcoded
        if ~exist('exp_name','var')
            % Ben data
            %exp_name='2015-03-03_AF03-light_awake';
            %exp_name='2015-03-04_AF11';
            %exp_name='2015-03-05_AF03';
            
            %exp_name='2015-03-05_AF11'; % example session 12Hz
            %exp_name='2015-03-05_AF11_compare';
            %exp_name='2015-03-05_AF11/resaved'; % recovered scim bitcodes!
            
            %exp_name='2015-04-07_AF11'; % example session 3Hz
            %exp_name='2015-04-07_AF11_compare';
            %exp_name='2015-04-10_AF11';
            
            %exp_name='2015-04-10_AF11_exp';
            %exp_name='2015-04-15_AF11'; % retinomapping
            %exp_name='2015-04-16_AF11';
            
            %exp_name='2015-07-17_AG02_awake';
            %exp_name='2015-07-20_AG02';
            %exp_name='2015-07-21_AG02';
            %exp_name='2015-07-22_AG02';
            
            
            
            %% 1024x300 sessions
            %exp_name='2015-08-03_AH02_init'; % no scim bitcodes...
            
            %exp_name='2015-08-03_AH02_init/resaved'; % recovered scim bitcodes by resaving in ImageJ
            %exp_name='2015-08-06_AH02';
            %exp_name='2015-08-10_AH02/resaved';
            %exp_name='2015-08-19_AH02'; % need ROIs
            
            
            %exp_name='2015-08-07_AH03'; % still downloading tif files
            %exp_name='2015-08-07_AH03/FOV01'; % new structure on NAS/Volume1/2photon/2p-data/raw
            
            
            %exp_name='2015-08-10_AH03';
            %exp_name='2015-08-14_AH03'; % session bitcodes needs cleaning up, did we do this? nopes
            %exp_name='2015-08-19_AH03';
            %exp_name='2015-08-21_AH03'; % funky eyedrift up
            
            %exp_name='2015-08-14_AH05'; % 
            %exp_name='2015-08-20_AH05'; % 
            %exp_name='2015-09-01_AH05'; % fix bitcodes: lsb was always off
            %exp_name='2015-09-02_AH05'; % no bitcodes... constant 15
            
            
            %exp_name='2015-08-18_AH06'; % fix bit Codes
            %exp_name='2015-08-26_AH06'; % red session mRuby2
            
            %exp_name='2015-09-01_AJ01';
            
            %exp_name='082715 SR101 vessels imaging'; % Ahbi vessels
            
            
            
            % Julie's data
            %exp_name='20150502_jrat3/Session02';
            %exp_name='20150502_jrat3/Session03-05';
            %exp_name='20150502_jrat3/Session06-10';
            %exp_name='20150502_jrat3/Session11-16';
            %exp_name='20150502_jrat3/Session17-22';
            %exp_name='20150502_jrat3/Session23-25';
            %exp_name='20150502_JR0013';
            
            
            %exp_name='20160321_JR005W'; % bitcode issues
            
            % Vincent lab data
            % ETL
            %exp_name='151218_KS154_etl_2P_KS/run03_ori8_reversed/plane01'; % werkt
            %X exp_name='151218_KS154_etl_2P_KS/run03_ori8_reversed_plane02'; % ugly
            %X exp_name='151218_KS154_etl_2P_KS/run03_ori8_reversed2_plane01'; % looks a lot like plane02 reversed
            %exp_name='151218_KS154_etl_2P_KS/run03_ori8_reversed2/plane02'; % looks like plane 01 reversed
            %exp_name='151218_KS154_etl_2P_KS/run03_movingbars_cardinal/plane01'; % plane 01!
            %exp_name='151218_KS154_etl_2P_KS/run03_movingbars_cardinal/plane02'; % plane 01 again... no switches?
            
            % cell
            %exp_name='150119_XH039_2P_XH\sf1_8dir_full';
            
            %exp_name='150206_MS049_2P_MS\run01_checkers2048'; % memory errors
            %exp_name='150122_KS127_2P_KS\run03_sf_tf_V1'; % lot's of processes
            %exp_name='150123_KS127_2P_KS\run02_retinotopy_V1';
            
            %exp_name='141221_DM044_2P_DM\run02_blank';
        end
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
%cd(root_folder)


% suppress output so we can run in -nojvm mode
animal_ID=exp_name(12:15);
%session_date=datevec(exp_name(1:11));
plot_it=0;
save_it=1;
use_GPU=0;
use_custom_ROI_defs=1; % Relevant for step 3 (GUI)

exp_name_txt=strrep(exp_name,'_',' ');
exp_name_txt=strrep(exp_name_txt,'\',': ');


% display current working folder
fprintf('Current working folder is : %s\n',data_folder)


%%% Define users
if use_custom_ROI_defs==1
    switch user_name
        case 'benvermaercke'
            ROI_definition_nr=2;
        case 'ben'
            ROI_definition_nr=2; % on the server | bkrunch
        case 'ben-pc\lbp'
            ROI_definition_nr=2; % KUL workstation
            
        case 'juliana'
            ROI_definition_nr=3;
            
        case 'julianarhee' % on the server
            ROI_definition_nr=3;
        case 'coxlab-2p\labuser'
            ROI_definition_nr=2;
            
            %case '' % template for other named users
            %    handles.ROI_definition_nr=4;
            
        otherwise
            ROI_definition_nr=1;
    end
else % store in
    ROI_definition_nr=1;
end