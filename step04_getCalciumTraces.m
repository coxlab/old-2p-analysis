clear all
clc

% 2DO
% how to handle shift in the middle of a file: motion correction in step 2a
% should solve that:
% BV20150410: have some version of motion correction, not completely robust
% yet
%
% BV20150409: add options to saved file + info about axis labels and such

header_script

session_vector=[];
switch 2
    case 1 % load single file recorded at same FOV, same stimuli
        cd(data_root)
        [file_name, pathname]=uigetfile('.mat','Pick data file');
        loadName=fullfile(pathname,file_name);
        %%
        load(loadName,'session_data')
        nSessions=1;
        count=1;
        valid_sessions(count)=session_data.data(1);
        valid_session_names{count}=loadName;
        valid_session_names_short{count}=file_name;
        session_vector=[6];
    case 2 % enitre folder and choose session
        data_folder=uigetdir(data_root);
        loadName=fullfile(data_folder,'data_analysis','session_overview.mat');
        load(loadName,'data_sessions')
        nSessions=length(data_sessions);
        
        %%% select only files that contain valid ROI_definitions
        count=1;
        valid_sessions=[];
        for iSess=1:nSessions
            [folder,file_name]=fileparts(data_sessions(iSess).file_name);
            loadName=fullfile(folder,'data_analysis',[file_name '.mat']);
            load(loadName,'session_data');
            if length(session_data.ROI_definitions)>1
                valid_sessions(count)=data_sessions(iSess).data(1);
                valid_session_names{count}=loadName;
                valid_session_names_short{count}=file_name;
                count=count+1;
            end
        end
        session_vector=valid_sessions;
        nSessions=length(session_vector);
    case 3 % specified folder and specified sessions
        data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-05_AF11/data_analysis';
        loadName_format='20150305_AF11_%03d.mat';
        session_vector=[1:4 8:12 15 17 18];
    case 4 % specified folder and specified sessions
        data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-15_AF11/data_analysis';
        loadName_format='20150305_AF11_%03d.mat';
        session_vector=[6];
        nSessions=length(session_vector);
end
%%
session_vector

%% construct a set of options that select different analysis steps to run on the data
% e.g. do motion correction, take out slow drifts, add f0 back in, calc z-score...
% add switches to select a scenario that replicates approaches from the
% literature, so we can easily compare
%options.motion_correction.use               =  0; % 0:not 1:TurboReg should be done before ROI definition

options.neuropil_subtraction.use            =  1;
options.neuropil_subtraction.factor         =  0.7;
options.neuropil_subtraction.TH_no_data     =  -4; % z-score

options.drift_correction.use                =  1;
options.drift_correction.averaging_interval =  15; % seconds
options.drift_correction.prctile            =  8;
options.drift_correction.add_mean           =  0;

options.calc_delta_f.method                 =  6;

options.save_data                           =  1;
options.plot_traces                         =  1;


%%
for iSess=1:nSessions
    session_nr=session_vector(iSess);
    loadName=valid_session_names{iSess};
    short_name=valid_session_names_short{iSess};
    short_name=strrep(short_name,'_',' ');
    load(loadName)
    loadName
    
    %%% Get all data from tif file
    tic
    frame_rate=session_data.data(5);
    file_name=session_data.file_name;
    info=imfinfo(file_name);
    nFrames=length(info);
    rows=session_data.data(4);
    cols=session_data.data(3);
    
    frames=zeros(rows-1,cols,nFrames);
    %mean_lum=zeros(nFrames,1);
    for iFrame=1:nFrames
        frame=double(imread(session_data.file_name,iFrame,'info',info));
        frames(:,:,iFrame)=double(frame(1:end-1,:)); % no flyback line double
        %mean_lum(iFrame)=mean(mean(frame(1:end-1,:)));
    end
    fprintf('Loading frames took: %3.1fs\n',toc);
    
    %% Extract ROI coordinates and construct time-series for each ROI
    % BV20150409: use motion correction values here to shift ROIs with the
    % movie.
    if isfield(session_data,'motion_correction')
        motion_correction=session_data.motion_correction;
        apply_motion_correction=1;
    else
        apply_motion_correction=0;
    end
    ROIs=session_data.ROI_definitions;
    nROI=length(ROIs);
    ROI_vector=cat(1,ROIs.ROI_nr);
    
    tic
    for iFrame=1:nFrames
        if apply_motion_correction==1
            offset=motion_correction.shift_matrix(iFrame,4:5);
        end
        
        for iROI=1:nROI
            if apply_motion_correction==1
                ROI_box=frames(offset(1)+ROIs(iROI).ROI_rect(2):offset(1)+ROIs(iROI).ROI_rect(4),offset(2)+ROIs(iROI).ROI_rect(1):offset(2)+ROIs(iROI).ROI_rect(3),iFrame); % do motion correction on ROIs
            else
                ROI_box=frames(ROIs(iROI).ROI_rect(2):ROIs(iROI).ROI_rect(4),ROIs(iROI).ROI_rect(1):ROIs(iROI).ROI_rect(3),iFrame); % do motion correction on ROIs
            end
            ROIs(iROI).timeseries_soma(iFrame)=mean(ROI_box(ROIs(iROI).mask_soma==1));
            ROIs(iROI).timeseries_neuropil(iFrame)=mean(ROI_box(ROIs(iROI).mask_neuropil==1));
        end
    end
    fprintf('Extracting %d ROIs data took: %3.1fs\n',[nROI toc]);
    
    %%% Construct activity matrix by processing previously extracted time-series
    activity_matrix=zeros(nROI,nFrames);
    tic
    
    
    %% Calculate number of dark frames (laser not on) at start of experiment
    if isfield(session_data,'blank_frames')
        blank_frames=session_data.blank_frames;
    else
        blank_frames=zscore(mean_lum)<options.neuropil_subtraction.TH_no_data;
        blank_frames(1:find(blank_frames>0,1,'last')+1)=true;
    end
    %sum(blank_frames)
    
    %%
    for iROI=1:nROI % for each ROI
        if options.neuropil_subtraction.use==0
            F=ROIs(iROI).timeseries_soma; % use raw soma signal
        else
            %%% Get extracted timeseries for soma and neuropil masks
            F_raw=ROIs(iROI).timeseries_soma;
            F_neuropil=ROIs(iROI).timeseries_neuropil;
            
            if any(blank_frames)
                F_raw(blank_frames)=mean(F_raw);
                F_neuropil(blank_frames)=mean(F_neuropil);
            end
            
            %%% subtract neuropil signals
            F=F_raw-F_neuropil*options.neuropil_subtraction.factor; % values .5 (Feinberg) and .7 (Kerlin) have been reported: .7 is used for 16x Nikon
        end
        
        %%% take out slow drifts
        if options.drift_correction.use==0
            F_no_drift=F; % leave signal as is
        else
            % calc number of frames that correspond to specified averaging interval
            nFrames_avg=options.drift_correction.averaging_interval*round(frame_rate);
            
            % make sure this number even
            if rem(nFrames_avg,2)==1
                nFrames_avg=nFrames_avg+1;
            end
            
            % take nFrames_avg around each sample and get 8th prctile value
            % of these values.
            prctile_vector=zeros(size(F));
            for iSample=1:nFrames
                if between(iSample,[nFrames_avg/2+1 nFrames-nFrames_avg/2])
                    hist_values=F(iSample-nFrames_avg/2+1:iSample+nFrames_avg/2);
                elseif iSample<nFrames_avg % handles start of trace
                    hist_values=F(1:iSample+nFrames_avg/2);
                elseif iSample>nFrames-nFrames_avg % handles end of trace
                    hist_values=F(iSample-nFrames_avg/2+1:end);
                end
                prctile_vector(iSample)=prctile(hist_values,options.drift_correction.prctile);
            end
            
            % subtract prctile vector from signal
            F_no_drift=F-prctile_vector;
            
            % optional, add average prctile value back in to restore mean
            % of trace and allow to check for fold changes (Feinberg 2015)
            % pointless when doing mean subtraction after this step
            if options.drift_correction.add_mean==1
                F_no_drift=F_no_drift+mean(prctile_vector);
            end
        end
        
        %%% Calculate delta F over F
        switch options.calc_delta_f.method
            case 0
                delta_F_no_drift=F_no_drift; % leave signal as is
                y_label='F';
                fixed_y_scale=30000;
            case 1 % naive way
                delta_F_no_drift=(F_no_drift-mean(F_no_drift))/mean(F_no_drift);
                y_label='\DeltaF/F';
                fixed_y_scale=30;
            case 2 % Feinberg
                delta_F_no_drift=F_no_drift/mean(F_no_drift)-1;
                y_label='\DeltaF/F';
                fixed_y_scale=30;
            case 3 % z-score naive
                delta_F_no_drift=(F_no_drift-mean(F_no_drift))/std(F_no_drift);
                y_label='Z(F)';
                fixed_y_scale=30;
            case 4 % minimal std for z-score
                disp('Under construction')
                delta_F_no_drift=(F_no_drift-median(F_no_drift));
                
                % estimate sigma, minimum from 10s sliding window sigma's
                averaging_interval=10; % seconds
                nFrames_avg=averaging_interval*round(frame_rate);
                if rem(nFrames_avg,2)==1
                    nFrames_avg=nFrames_avg+1;
                end
                
                sigma_vector=zeros(size(delta_F_no_drift));
                for iSample=1:nFrames
                    if between(iSample,[nFrames_avg/2+1 nFrames-nFrames_avg/2])
                        sigma_values=delta_F_no_drift(iSample-nFrames_avg/2+1:iSample+nFrames_avg/2);
                    elseif iSample<nFrames_avg
                        sigma_values=delta_F_no_drift(1:iSample+nFrames_avg/2);
                    elseif iSample>nFrames-nFrames_avg
                        sigma_values=delta_F_no_drift(iSample-nFrames_avg/2+1:end);
                    else
                        die
                    end
                    sigma_vector(iSample)=std(sigma_values);
                end
                sigma_est=min(sigma_vector);
                delta_F_no_drift=delta_F_no_drift/sigma_est;
                y_label='Z(F)';
                fixed_y_scale=30;
            case 5 % z-score based on e-phys analysis
                % find regions without 'spike' and use rest for std and mean calculation
                F_temp=F_no_drift-median(F_no_drift);
                noiseSTD=median(abs(F_temp))/.6745;
                min_heigth_factor=4;
                TH=noiseSTD*min_heigth_factor;
                selection_vector=F_temp<TH;
                mu=mean(F_no_drift(selection_vector));
                sigma=std(F_no_drift(selection_vector));
                
                delta_F_no_drift=(F_no_drift-mu)/sigma;
                if 0
                    %%
                    T=1:length(F_temp);
                    plot(T,F_temp)
                    hold on
                    plot([0 length(F_temp)],[TH TH])
                    plot(T(selection_vector==1),F_temp(selection_vector==1),'r')
                    hold off
                end
                y_label='Z(F)';
                fixed_y_scale=30;
            case 6 % deltaF over F using non-modulation parts
                % use rolling average to select no-modulation parts
                F_no_drift_smooth=smooth(F_no_drift,round(frame_rate*4));
                
                TH=prctile(F_no_drift_smooth,2)+512; % optimize!!!
                selection_vector=F_no_drift_smooth<TH;
                mu=mean(F_no_drift(selection_vector));
                sigma=std(F_no_drift(selection_vector));
                
                delta_F_no_drift=(F_no_drift-mu)/mu;
                y_label='\DeltaF/F';
                fixed_y_scale=30;
                
                if 0
                    %%
                    F_temp=F_no_drift_smooth;
                    T=1:length(F_temp);
                    plot(T,F_temp)
                    hold on
                    %plot(T,F_no_drift)
                    plot([0 length(F_temp)],[TH TH])
                    plot(T(selection_vector==1),F_temp(selection_vector==1),'r.','MarkerSize',7)
                    hold off
                end
        end
        options.y_label=y_label;
        options.fixed_y_scale=fixed_y_scale;
        
        %%% Replace blank frames at start of movie by average of final trace
        %delta_F_no_drift(sel_no_data)=mean(delta_F_no_drift);
        delta_F_no_drift(blank_frames)=0;
        
        if 0
            %%
            plot(delta_F_no_drift)
            hold on
            plot([0 length(delta_F_no_drift)],[0 0],'r')
            hold off
            ylabel(y_scale)
            axis([0 length(delta_F_no_drift) -fixed_y_scale*.1 fixed_y_scale])
            box off
        end
        
        %%% Stored processed trace to activity matrix
        if any(isnan(delta_F_no_drift))
            die
        end
        activity_matrix(iROI,:)=delta_F_no_drift;
        
    end
    fprintf('Constructing acitivity matrix took: %3.1fs\n',toc);
    
    if options.save_data==1
        %% Save activity and stim_matrix to file
        session_data.activity_matrix=activity_matrix';
        session_data.options=options;
        %session_data.mean_lum=mean_lum;
        save(loadName,'session_data')
    end
    
    nROI
    if options.plot_traces==1
        %% plot all traces
        figure(1)
        clf
        
        cols=ceil(sqrt(nROI));
        rows=ceil((nROI+1)/cols);
        stim_data=session_data.stimulus_matrix_ext;
        stim_vector=stim_data(:,4)>0;
        for iROI=1:nROI
            subplot(rows,cols,iROI)
            ROI_nr=ROIs(iROI).ROI_nr;
            
            F=activity_matrix(iROI,:);
            max_val=max([fixed_y_scale max(F)]);
            prop=session_data.ROI_definitions(iROI).ellipse_properties;
            info=sprintf('Radius: %3.2f',prop.long_axis);
            
            bar(spherify(stim_vector,2)*max_val,'barWidth',1,'FaceColor',[1 1 1]*.9,'EdgeColor',[1 1 1]*.9)
            hold on
            plot(F,'color',[0 1 0]*.7)
            text(nFrames/10,-max_val*.1,info)
            hold off
            
            ylabel(options.y_label)
            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            axis([1 nFrames -max_val*.2 max_val])
            
            title({short_name,sprintf('Session %02d-Cell #%03d',[session_nr ROI_nr])})
            drawnow
        end
        
        %%
        CC=corr(activity_matrix');
        CC_between=CC(eye(nROI)==0);
        [r,c]=find(CC==max(CC_between));
        [ROI_vector(r(1)) ROI_vector(c(1)) max(CC_between)]
        subplot(rows,cols,iROI+1)
        imagesc(CC)
        set(gca,'clim',[-1 1])
        colorbar
        set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
        axis square
        
        %%
        frame_selection=2:min([nFrames size(session_data.stimulus_matrix_ext,1)]);
        
        stim_data=session_data.stimulus_matrix_ext(frame_selection,:);
        stim_frames=stim_data(:,4)>0;
        
        figure(2)
        imagesc(activity_matrix(:,frame_selection))
        hold on
        plot(spherify(-stim_vector,2)+nROI-3,'r')
        hold off
        colormap(green)
    end
end

disp('All done')