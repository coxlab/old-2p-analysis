clear all
clc


header_script
%data_root='/Users/benvermaercke/Dropbox (coxlab)/2p-data/';
%exp_folder='2015-07-17_AG02_awake';
%data_folder=fullfile(data_root,exp_folder);
files=scandir(data_folder,'tif');
nFiles=length(files);


for iFile=7%:nFiles
    file_name=fullfile(data_folder,files(iFile).name);
    if exist(file_name,'file')==2
        fprintf('Pre-processing file %s...\n',file_name)
        save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
        save_name=strrep(save_name,'tif','mat');
        
        if exist(save_name,'file')==0
            %%% Create file based on class file
            session_data=imaging_dataset(file_name);
            session_data.save_data()
        else
            load(save_name,'session_data')
        end
        
        %%% Extract info from movie
        session_data.get_mov_info()
        session_data.get_scim_data()
        session_data.read_flyback()
        session_data.get_FOV_info(.85)
        
        %%% Detect invalid frames
        session_data.find_blank_frames() % due to laser power not being on
        
        if session_data.is_static_FOV()==0 % if recording at 1 position
            fprintf('Session does not appear to be a static FOV...\n')
        else            
            %%% Extract bitcode information from two sources
            session_data.get_scim_bitCodes()
            session_data.get_MWorks_bitCodes()
            session_data.find_offset()
            fprintf('Offset was determined to be %d events...\n',session_data.bitCodes.offset)
            
            %%% Save results so far
            session_data.save_data()
            
            %%% Motion Correction
            session_data.find_reference_image()
            if 0
                %%
                imshow(calc_gamma(session_data.motion_correction.reference_image.im,.5),[])
                colormap(green)
            end
            session_data.do_motion_correction()
            session_data.find_motion_frames(2)
            if 0
                %%
                plot(session_data.motion_correction.shift_matrix(:,2))
                hold on
                bar(double(session_data.motion_correction.ignore_frames))
            end
            
            if 0
                %%
                session_data.visualize_motion_correction(1)
            end
            
            %%% Calculate z-projections
            session_data.do_calc_MIPs()
            if 1
                %%
                session_data.imshow(session_data.MIP_std.data)
            end
            session_data.save_data()
            
            %session_data.do_motion_detection()
            % plot(session_data.motion_proxy)
        end
    end
end

%%
%%% After all preprocessing, compile session overview file so we can run
%%% manual ROI definition

% Run step 03
% Rest of pipeline is sort of same
for iFile=6%:nFiles
    save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
    save_name=strrep(save_name,'tif','mat');
    load(save_name,'session_data') % reload after step03, probably needs to be separate script
    %if session_data.is_static_FOV()&&~isempty(fieldnames(session_data.ROI_definitions))
    %ROI_definition_nr=2;
    if length(session_data.ROI_definitions)==ROI_definition_nr&&~isempty(session_data.ROI_definitions(ROI_definition_nr).ROI(1).ROI_nr)
        %%% Extract activity traces
        %session_data.reset_trace_matrix() % allows to recalculate the traces
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
        
        %%
        session_data.combine_act_stim(4,6)
        
    else
        step03_ROI_GUI()
    end
end




if 0
    %% get stack
    save_name=fullfile(data_folder,'data_analysis',files(5).name);
    save_name=strrep(save_name,'tif','mat');
    load(save_name,'session_data')
    stack=session_data;
    
    for iFile=1%:nFiles
        save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
        save_name=strrep(save_name,'tif','mat');
        load(save_name,'session_data')
        if session_data.is_static_FOV()
            session_data.find_FOV_in_stack(stack);
        else % is stack? find other field in this stack if same xy location, give z
            xyz=cat(1,session_data.frame_info.xyz_submicron);
            Z=xyz(:,3);
            plot(Z)
            
        end
    end
    
end

