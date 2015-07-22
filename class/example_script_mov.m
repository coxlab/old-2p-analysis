clear all
clc

file_name='/Users/benvermaercke/Dropbox (coxlab)/GRat06/150716/rat7.tif';

session_data=imaging_dataset(file_name);
session_data.save_data()

%%
session_data.get_mov_info()
session_data.mov_info.frame_rate=30;
session_data.find_blank_frames()
%%
session_data.find_reference_image()

%%
session_data.do_motion_correction()
session_data.find_motion_frames(2)
%%
session_data.do_calc_MIPs()
session_data.save_data()

%% do step 03 here

%%
session_data.load_data()
ROI_definition_nr=2

session_data.Activity_traces.extraction_options.neuropil_subtraction.use=0;
session_data.Activity_traces.extraction_options.drift_correction.averaging_interval=2;
session_data.do_trace_extraction(ROI_definition_nr)

figure(1)
plot(session_data.Activity_traces.activity_matrix)
