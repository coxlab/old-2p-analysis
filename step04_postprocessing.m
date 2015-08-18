clear all
clc

header_script


files=scandir(data_folder,'tif');
nFiles=length(files);
if nFiles==0
    % In case we moved the raw files to a subfolder to allow those to be
    % selectively unsynced
    data_folder=fullfile(data_folder,'tif_files');
    files=scandir(data_folder,'tif');
    nFiles=length(files);
end
%%
%%% After all preprocessing, compile session overview file so we can run
%%% manual ROI definition

% Run step 03
% Rest of pipeline is sort of same

%ROI_definition_nr=1; % use auto ROIs

for iFile=1:nFiles
    save_name=fullfile(data_folder,'data_analysis',files(iFile).name)
    save_name=strrep(save_name,'tif','mat');
    load(save_name,'session_data') % reload after step03, probably needs to be separate script
    %if session_data.is_static_FOV()&&~isempty(fieldnames(session_data.ROI_definitions))
    %ROI_definition_nr=2;
    
    if session_data.is_static_FOV==1
        session_data.ROI_definition_nr=ROI_definition_nr;
        if length(session_data.ROI_definitions)>=ROI_definition_nr&&~isempty(session_data.ROI_definitions(ROI_definition_nr).ROI(1).ROI_nr)
            
            
            session_data.ROI_definition_nr=ROI_definition_nr;
            
            %%% Extract activity traces
            %session_data.reset_trace_matrix() % allows to recalculate the traces
            session_data.do_trace_extraction()
            session_data.save_data()
            %session_data.plot_traces()
            
            %% Extract stimulus relevant information
            %session_data.bitCodes.MWorks_bitCodes=[];
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
                rethrow(lasterror)
            end
           
        else
            %step03_ROI_GUI()
            die
        end
    end
end