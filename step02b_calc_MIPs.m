clear all
clc

% BV20150416: apply new idea for motion correction:
% - use frame to frame shift to detect homogeneous blocks (=noisy)
% - average over blocks and detect big shifts (=reliable)
% - compute offset between block average and single frames within block (=better reliability then step 01) 
%
% BV20150420: prepend median filter 1-1-3 to MIP calculations, won't affect
% raw file, so don't care about temporal smoothing this causes.

header_script

debugging=0;
MIP_options.do_medfilt3=1;
MIP_options.medfilt_values=[1 1 3]; % default 1-1-3, could try 1-1-5
MIP_options.apply_motion_correction=1; % try only if shift_values are present

%%% Use uigetdir to get data folder
%cd(data_root)
%data_folder=uigetdir(data_root);

loadName=fullfile(data_folder,'data_analysis','session_overview.mat');
load(loadName,'data_sessions')
nSessions=length(data_sessions);

%%
t0=clock;
for iSess=1:nSessions
    [folder,file_name]=fileparts(data_sessions(iSess).file_name);
    loadName=fullfile(folder,'data_analysis',[file_name '.mat']);
    load(loadName,'session_data');
    
    info=imfinfo(session_data.file_name);
    nFrames=session_data.data(2);
    rows=session_data.data(4);
    cols=session_data.data(3);
    
    frames=zeros(rows-1,cols,nFrames);
    if isfield(session_data,'motion_correction')
        motion_correction=session_data.motion_correction;
        MIP_options.apply_motion_correction=1;
    else
        MIP_options.apply_motion_correction=0;
    end
    
    fprintf('Loading frames...')
    tic
    for iFrame=1:nFrames
        frame=double(imread(session_data.file_name,iFrame,'info',info));
        frame=double(frame(1:end-1,:)); % no flyback line
        
        if MIP_options.apply_motion_correction==1
            %%% BV20150409: use motion correction here to offset frames
            % will create boundary issues
            offset=-motion_correction.shift_matrix(iFrame,4:5);
            if any(offset~=0)
                %frame=offsetIm(frame,offset(1),offset(2),0);
                frame=offsetIm(frame,offset(1),offset(2),mean(frame(:)));
            else
                % as long as offsets are both zero, use uncorrected version
            end
        end
        frames(:,:,iFrame)=frame;
    end
    fprintf('took %3.2f seconds.\n',toc)
    if MIP_options.do_medfilt3==1&&debugging==0
        %% BV20150420: added 3D median filt, will make images look prettier, especially max
        tic
        frames=medfilt3(frames,MIP_options.medfilt_values);
        fprintf('Medfilt3(A,[%d %d %d]) took %3.2f seconds.\n',[MIP_options.medfilt_values toc])
    end    
    
    %% calc projection images
    tic
    fprintf('Computing MIPs...')
    session_data.MIP_avg.data=mean(frames,3);
    session_data.MIP_avg.gamma_val=1;
    session_data.MIP_max.data=max(frames,[],3);
    session_data.MIP_max.gamma_val=.6;
    session_data.MIP_std.data=std(frames,[],3);
    session_data.MIP_std.gamma_val=.8;
    if debugging==0 % only if we are not debugging, takes a long time to run
        session_data.MIP_cc_local.data=CrossCorrImage(frames);
        session_data.MIP_cc_local.gamma_val=.6;
    end
    fprintf('took %3.2f seconds.\n',toc)
    
    if 0
        %% cross-correlation of most active pixels : under construction
        tic
        session_data.MIP_cc_global.data=CrossCorrImage_global(HandleObject(frames),session_data.MIP_std.data);
        session_data.MIP_cc_global.gamma_val=1;
        toc
    end
    
    if debugging==1
        %%
        subplot(221)
        imshow(calc_gamma(session_data.MIP_avg.data,session_data.MIP_avg.gamma_val),[])
        subplot(222)
        imshow(calc_gamma(session_data.MIP_max.data,session_data.MIP_max.gamma_val),[])
        subplot(223)
        imshow(calc_gamma(session_data.MIP_std.data,session_data.MIP_std.gamma_val),[])
        colormap(green)
    end
    %session_data.frames=frames;
    
    if debugging==0
        %%
        session_data.MIP_options=MIP_options;
        [save_folder, save_name]=fileparts(session_data.file_name);
        saveName=fullfile(save_folder,'data_analysis',[save_name '.mat'])
        save(saveName,'session_data')
        progress(iSess,nSessions,t0)
    end
end