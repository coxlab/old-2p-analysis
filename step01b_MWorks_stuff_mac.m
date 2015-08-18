clear all
clc

header_script

files=scandir(fullfile(data_folder,'data_analysis'),'mat');
nFiles=length(files);

for iFile=1:nFiles
    load_name=fullfile(fullfile(data_folder,'data_analysis'),files(iFile).name);
    load(load_name,'session_data')
    try
        session_data.rebase(data_root)
        ROI_definition_nr
        session_data.ROI_definition_nr=ROI_definition_nr;
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