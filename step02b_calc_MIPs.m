clear all
clc

header_script

%%% Use uigetdir for general use
cd(data_root)
data_folder=uigetdir(data_root);

loadName=fullfile(data_folder,'data_analysis','session_overview.mat');
load(loadName,'data_sessions')
nSessions=length(data_sessions);

%%
t0=clock;
for iSess=2:nSessions
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
        apply_motion_correction=1;
    else
        apply_motion_correction=0;
    end
    tic
    for iFrame=1:nFrames
        frame=double(imread(session_data.file_name,iFrame,'info',info));
        frame=double(frame(1:end-1,:)); % no flyback line
        
        if apply_motion_correction==1
            %%% BV20150409: use motion correction here to offset frames
            % will create boundary issues
            offset=-motion_correction.shift_matrix(iFrame,4:5);
            if any(offset~=0)
                %frame=offsetIm(frame,offset(1),offset(2),0);
                frame=offsetIm(frame,offset(1),offset(2),mean(frame(:)));
                die
            else
                % as long as offsets are both zero, use uncorrected version
            end
        end
        frames(:,:,iFrame)=frame;
    end
    fprintf('Loading frames took %3.2f seconds.\n',toc)
    
    %% calc projection images
    session_data.MIP_avg.data=mean(frames,3);
    session_data.MIP_avg.gamma_val=1;
    session_data.MIP_max.data=max(frames,[],3);
    session_data.MIP_max.gamma_val=1;
    session_data.MIP_std.data=std(frames,[],3);
    session_data.MIP_std.gamma_val=1;
    if 1 % only if we are not debugging, takes a long time to run
        session_data.MIP_cc_local.data=CrossCorrImage(frames);
        session_data.MIP_cc_local.gamma_val=1;
    end
    
    if 0
        %% cross-correlation of most active pixels : under construction
        tic
        session_data.MIP_cc_global.data=CrossCorrImage_global(HandleObject(frames),session_data.MIP_std.data);
        session_data.MIP_cc_global.gamma_val=1;
        toc
    end
    
    if 0
        %%
        subplot(221)
        imshow(calc_gamma(session_data.MIP_avg.data,session_data.MIP_avg.gamma_val),[])
        subplot(222)
        imshow(calc_gamma(session_data.MIP_max.data,session_data.MIP_max.gamma_val),[])
        subplot(223)
        imshow(calc_gamma(session_data.MIP_std.data,session_data.MIP_std.gamma_val),[])
        %subplot(224)
        %imshow(calc_gamma(session_data.MIP_cc_local.data,session_data.MIP_cc_local.gamma_val),[])
        colormap(green)
    end
    %session_data.frames=frames;
    
    if 1
        %%
        [save_folder, save_name]=fileparts(session_data.file_name);
        saveName=fullfile(save_folder,'data_analysis',[save_name '.mat'])
        save(saveName,'session_data')
        progress(iSess,nSessions,t0)
    end
end