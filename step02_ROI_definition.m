clear all
clc

header_script
%data_root='/Users/benvermaercke/Dropbox (coxlab)/2p-data/';
%exp_folder='2015-07-17_AG02_awake';
%data_folder=fullfile(data_root,exp_folder);
files=scandir(data_folder,'tif');
nFiles=length(files);
%%

for iFile=1:nFiles
    file_name=fullfile(data_folder,files(iFile).name);
    if exist(file_name,'file')==2
        fprintf('Finding ROIs for file %s...\n',file_name)
        save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
        save_name=strrep(save_name,'tif','mat');
        
        if exist(save_name,'file')==0
            error('Run step 01 first...')
        else
            load(save_name,'session_data')
            
            
            if session_data.is_static_FOV==1
                %%% Make sure filenames are relative to data folder on this machine
                session_data.rebase(data_root)
                
                %%% Calc ROIs based on adaptive thresholding
                %session_data.reset_ROIs()
                session_data.find_ROIs()
                
                %%% Save to ROI_definitions, first field
                session_data.save_data()
                if 0
                    %%
                    session_data.ROI_definitions(1).ROI;
                    session_data.plot_ROIs()
                end
            end
        end
    end
end