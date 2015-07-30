%%
%%% After all preprocessing, compile session overview file so we can run
%%% manual ROI definition

% Run step 03
% Rest of pipeline is sort of same
for iFile=1:nFiles
    save_name=fullfile(data_folder,'data_analysis',files(iFile).name)
    save_name=strrep(save_name,'tif','mat');
    load(save_name,'session_data') % reload after step03, probably needs to be separate script
    %if session_data.is_static_FOV()&&~isempty(fieldnames(session_data.ROI_definitions))
    %ROI_definition_nr=2;
    if length(session_data.ROI_definitions)==ROI_definition_nr&&~isempty(session_data.ROI_definitions(ROI_definition_nr).ROI(1).ROI_nr)
        %%% Extract activity traces
        session_data.reset_trace_matrix() % allows to recalculate the traces
        session_data.do_trace_extraction(ROI_definition_nr)
        session_data.save_data()
        %session_data.plot_traces()
        
        %% Extract stimulus relevant information
        %session_data.bitCodes.MWorks_bitCodes=[];
        session_data.get_MWorks_bitCodes()
        session_data.get_exp_type()
        session_data.get_MWorks_stimulus_info()
        session_data.create_stim_matrix()
        %session_data.Experiment_info.stim_matrix
        session_data.save_data()
        %%
        %session_data.combine_act_stim(1,6)
        
    else
        step03_ROI_GUI()
        die
    end
end