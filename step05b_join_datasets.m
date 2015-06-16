clear all
clc

%%% BV20150507: init
%%% This script will take data from the same FOV in the same session (same
%%% day) and concatenate those into 1 dataset.
%%% This database file will form the basis for all subsequent analysis.
%%% The Stim and Resp matrices will be the main components and will allow
%%% for easy construction of RF by various methods, preceded by different
%%% kinds of deconvolution/mvregression.
%%% Use clusters defined in step 1b to join datafiles.
%%% The resulting database file will contain:
%%%     the joined stim and resp matrices
%%%     timescale
%%%     as well as MIPs and ROI locations for plotting
%%%     properties per traces such as #transients,
%%%     consider also adding stim averaged data

header_script
%save_it=1;

loadName=fullfile(data_folder,'data_analysis','session_overview.mat');
load(loadName,'data_sessions','FOV_matching')

nClusters=max(FOV_matching.clusters(:));
for iClust=1:nClusters
    [r,c]=find(FOV_matching.clusters==iClust);
    selected_sessions=unique(r);
    
    %%% Load data for selected sessions
    nSessions=length(selected_sessions);
    session_data_all=cell(nSessions,1);
    count=1;
    for iSess=1:nSessions
        [f,fn,ext]=fileparts(data_sessions(selected_sessions(iSess)).file_name);
        
        loadName=fullfile(data_folder,'data_analysis',[fn '.mat']);
        if exist(loadName,'file')
            load(loadName,'session_data')
            if isfield(session_data,'options')
                % remove old options field
                session_data=rmfield(session_data,'options');
            end
            if isfield(session_data,'activity_matrix')
                %session_data_all(count)=session_data;
                session_data_all{count}=session_data;
                count=count+1;
            end
        end
    end
    
    %%% Get all unique numbers for this FOV
    unique_cell_numbers=[];
    for iSess=1:nSessions
        session_data=session_data_all{iSess};
        ROIs=get_ROI_definitions(session_data,ROI_definition_nr);
        %R=session_data_all{iSess}.ROI_definitions;
        cell_numbers=cat(1,ROIs.ROI_nr);
        if isempty(unique_cell_numbers)
            unique_cell_numbers=cell_numbers;
        else
            unique_cell_numbers=intersect(unique_cell_numbers,cell_numbers);
        end
    end
    nROI=length(unique_cell_numbers);
    
    %%% Concatenate both stimulus and activity data
    STIM_ALL=[];
    RESP_ALL=[];
    RESP_ALL_oopsi=[];
    
    timescale_all=[];
    cur_time=0;
    for iSess=1:nSessions
        session_data=session_data_all{iSess};
        
        %%% Construct timescale
        nFrames=session_data.data(2);
        frame_rate=session_data.data(5);
        T=((1:nFrames)-1)/frame_rate;
        
        %%% Get data
        ROIs=get_ROI_definitions(session_data,ROI_definition_nr);
        cell_numbers=cat(1,ROIs.ROI_nr);
        cell_selection=ismember(cell_numbers,unique_cell_numbers);
        time_selection=1:size(session_data_all{iSess}.stimulus_matrix_ext,1);
        
        %%% clip off last trial if not ending in blank
        STIM=session_data.stimulus_matrix_ext(time_selection,:);
        switch session_data.expType
            case 1
                trial_vector=STIM(:,4);
            case 2
                trial_vector=STIM(:,5);
        end
        
        clipped_frames=false(nFrames,1);
        if trial_vector(end)>-1
            last_blank_frame=find(trial_vector==-1,1,'last');
            clipped_frames(last_blank_frame+1:end)=true;
        end
        
        % skip blank and motion-frenzy frames
        time_selection(session_data.blank_frames==1|session_data.motion_correction.ignore_frames==1|clipped_frames==1)=[];
        
        % skip unreliable motion frames
        %time_selection()=[];
        nFrames_selected=length(time_selection);
        
        % concat time scales
        timescale=T(time_selection)-min(T(time_selection));
        timescale_all=cat(1,timescale_all,cur_time+timescale');
        cur_time=cur_time+max(timescale);
        
        % concat stim and resp matrices
        STIM=session_data.stimulus_matrix_ext(time_selection,:);
        RESP=session_data.activity_matrix(time_selection,cell_selection);
        RESP_oopsi=session_data.spike_matrix(time_selection,cell_selection);
        
        STIM_ALL=cat(1,STIM_ALL,STIM);
        RESP_ALL=cat(1,RESP_ALL,RESP);
        RESP_ALL_oopsi=cat(1,RESP_ALL_oopsi,RESP_oopsi);
    end
    
    %%% Post hoc stuff
    % Get number of sign transients, some people turn the other parts of
    % the track to zero
    TH=5;
    nSpikes_vector=zeros(nROI,1);
    for iROI=1:nROI
        sel=RESP_ALL(:,iROI)>TH;
        start_points=find(diff(sel)==1);
        nSpikes_vector(iROI)=length(start_points);
    end
    
    [sorted,order]=sort(nSpikes_vector,'descend');
    cell_numbers_ranked=unique_cell_numbers(order);
    
    %%% Construct dataset
    dataset.cluster_nr=iClust;
    dataset.session_vector=selected_sessions;
    
    dataset.nFrames=size(RESP_ALL,1);
    dataset.timescale=timescale_all;
    
    dataset.stim_matrix=STIM_ALL;
    dataset.resp_matrix=RESP_ALL;
    dataset.resp_matrix_NND=RESP_ALL_oopsi;
    dataset.nSpikes_vector=nSpikes_vector;
    dataset.cell_numbers_ranked=cell_numbers_ranked;
    dataset.FOV_info=session_data.FOV_info;
    dataset.MIP=session_data.MIP_std;
    dataset.ROI_definitions=session_data.ROI_definitions;
    dataset.nROI=nROI;
    
    %%% Save data
    if save_it==1
        %%
        saveName=fullfile(data_folder,'data_analysis',sprintf('dataset_%03d.mat',iClust));
        if exist(saveName,'file')
            %dataset_update=dataset;
            %load(saveName,'dataset')
            % add what is new/changed
            %dataset.session_vector=dataset_update.session_vector;
            save(saveName,'dataset')            
        else
            save(saveName,'dataset')
        end
        fprintf('Saved dataset to %s.\n',saveName)
    end
end
