% Start from this script to create a full pipeline containing separate
% chapters: 
% - preprocessing - use cluster!
% - manual ROI definition - locally!
% - postprocessing/dataset creation
% - various kinds of data analysis : using different class!


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
        fprintf('Pre-processing file %s...\n',file_name)
        save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
        save_name=strrep(save_name,'tif','mat');
        
        if exist(save_name,'file')==0
            %%% Create file based on class file
            session_data=imaging_datasession(file_name);
            session_data.save_data()
        else
            load(save_name,'session_data')
        end
        
        %%% Make sure filenames are relative to data folder on this machine
        session_data.rebase(data_root)
        
        %%% Extract info from movie
        session_data.get_mov_info()
        session_data.get_scim_data()
        session_data.read_flyback()
        session_data.get_FOV_info(.85)
        
        %%% Detect invalid frames
        session_data.find_blank_frames() % due to laser power not being on
        
        session_data.save_data()
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
            %session_data.reset_motion_correction();
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
            if 0
                %%
                session_data.imshow(session_data.MIP_std.data)
            end
            session_data.save_data()
            
            %session_data.do_motion_detection()
            % plot(session_data.motion_proxy)
        end
    end
end
di

if 0
    %% 
    iFile=1
    file_name=fullfile(data_folder,files(iFile).name);
    if exist(file_name,'file')==2
        fprintf('Exporting file %s...\n',file_name)
        save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
        save_name=strrep(save_name,'tif','mat');
        load(save_name)
        session_data.export_movie('test_01.tif',session_data.get_frames([],1))
    end
    disp('Done!')
end

%%
%%% After all preprocessing, compile session overview file so we can run
%%% manual ROI definition

% Run step 03
% Rest of pipeline is sort of same
for iFile=1:nFiles
    save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
    save_name=strrep(save_name,'tif','mat');
    load(save_name,'session_data') % reload after step03, probably needs to be separate script    
    %if session_data.is_static_FOV()&&~isempty(fieldnames(session_data.ROI_definitions))
    %ROI_definition_nr=2;
    session_data.ROI_definition_nr=ROI_definition_nr;
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
        session_data.save_data()
        %%
        %session_data.combine_act_stim(1,6)
        
    else
        step03_ROI_GUI()
        die
    end
end

%% Join datasets

% check folder for joinable sessions: close in space (and time)
[clusters,cluster_vector,nClusters,files]=session_data.find_FOV_clusters(session_data.folder_info.save_folder);

% if close in space, but not in time, link files but call it a different
% timepoint, so you can look at evolution and/or stability

for iClust=1:nClusters
    cluster_nr=cluster_vector(iClust);
    session_vector=find(clusters==cluster_nr);
    nSessions=length(session_vector);
    
    %% Load all session into one compound object
    clear S
    for iSession=1:nSessions
        load_name=fullfile(session_data.folder_info.save_folder,files(iSession).name);
        load(load_name,'session_data')
        S(iSession)=session_data;
    end
    
    %% Check whether same number of ROIs are defined in all sessions
    ROI_counts=S.get_ROI_counts(ROI_definition_nr);
    if all(ROI_counts==ROI_counts(1))
        nROIs=mean(ROI_counts);
    else
        error('Number of ROIs has to be identical for all sessions');
    end
    
    
    %% Verify distances are closely matched over sessions
    distances=S.get_ROI_distances(ROI_definition_nr);
    if any(mean(distances)>10)
        error('Large shift in ROI centers detected...')
    else
        % all good: take the ROI_definitions of one of the sessions
        %ROI_definitions=S(1).ROI_definitions(ROI_definition_nr).ROI;
        
        % Then join them 
        % leave out last trial
        dataset=S.join_data_sessions(ROI_definition_nr);
        dataset.cluster_nr=iClust;
        dataset.session_vector=session_vector;
        dataset.nFrames=size(dataset.STIM,1);
    end
            
    %% Do something with dataset
    %dataset.plot_traces()
    dataset.RF_analysis(13)
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

