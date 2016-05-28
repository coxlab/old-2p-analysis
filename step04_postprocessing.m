clear all
clc

header_script

files=scandir(fullfile(data_folder,'data_analysis'),'mat');
nFiles=length(files);

% 
% files=scandir(data_folder,'tif');
% nFiles=length(files);
% tif_in_sub_folder=0;
% if nFiles==0
%     % In case we moved the raw files to a subfolder to allow those to be
%     % selectively unsynced
%     data_folder_new=fullfile(data_folder,'tif_files');
%     files=scandir(data_folder_new,'tif');
%     nFiles=length(files);
%     tif_in_sub_folder=1;
% end
%%
%%% After all preprocessing, compile session overview file so we can run
%%% manual ROI definition

% Run step 03
% Rest of pipeline is sort of same



%ROI_definition_nr=1; % use auto ROIs

for iFile=1:nFiles
    
    %save_name=fullfile(data_folder,'data_analysis',files(iFile).name)
    %save_name=strrep(save_name,'tif','mat');
    load_name=fullfile(data_folder,'data_analysis',files(iFile).name);
    load(load_name,'session_data') % reload after step03, probably needs to be separate script
    %if session_data.is_static_FOV()&&~isempty(fieldnames(session_data.ROI_definitions))
    %ROI_definition_nr=2;
    
    %%% Make sure paths are relative to our current Dropbox location
    session_data.rebase(data_root)
    if isdir(data_folder)&&~isempty(scandir(data_folder,'tif'))
        session_data.rebase_tif(data_folder)
    elseif isdir(fullfile(data_folder,'tif_files'))
        session_data.rebase_tif(fullfile(data_folder,'tif_files'))
    else
        error('Raw tif files not found, need those for ROI extraction...')
    end
        
    if session_data.is_static_FOV==1&&session_data.mov_info.nFrames>300
        session_data.ROI_definition_nr=ROI_definition_nr;        
        if length(session_data.ROI_definitions)>=ROI_definition_nr&&~isempty(session_data.ROI_definitions(ROI_definition_nr).ROI(1).ROI_nr)
            
            %session_data.ROI_definition_nr=ROI_definition_nr;
            
            
            %%% Create sparse ROI regions based on coords_MIP
            tic
            session_data.create_mask_from_ROI()
            session_data.clean_neuropil_shell()
            toc
            
            %%% Extract activity traces
            session_data.reset_trace_matrix() % allows to recalculate the traces
            session_data.Activity_traces.extraction_options.calc_delta_f_method=0;
            session_data.do_trace_extraction()
            %session_data.save_data()
            %session_data.plot_traces()
            
            %% Extract stimulus relevant information
            %session_data.bitCodes.MWorks_bitCodes=[];
            if ismac==1 % don't try this stuff on the server, use step1b first
                %%
                try
                    session_data.get_scim_bitCodes()
                    session_data.get_MWorks_bitCodes()
                    session_data.find_offset()
                    session_data.get_MWorks_bitCodes()
                    session_data.get_exp_type()
                    session_data.get_MWorks_stimulus_info()
                    session_data.create_stim_matrix()
                    %session_data.Experiment_info.stim_matrix
                    session_data.save_data()
                catch
                    E=lasterror;
                    disp(E.message)
                    disp('This file was not analysed properly')
                    %rethrow(lasterror)
                end
            end
            
        else
            error('No custom ROIs defined... Run Step03 first!!!')
            %step03_ROI_GUI()
            
        end
    end
end