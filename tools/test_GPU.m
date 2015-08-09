clear all
clc

header_script
load(fullfile(data_folder,'data_analysis','2015-08-07_AH03_001.mat'))

session_data.rebase(data_root)

if 0
    frames=double(session_data.get_frames());
    
    ref=mean(frames(:,:,20:30),3);
    nFrames=size(frames,3);
end

if 0
    %% Single frame case
    disp('1) Single frame case, conversion within function...')
    tic
    [CC_max,offset]=im_align(frames(:,:,10),ref,0);
    T1=toc;
    
    tic
    [CC_max,offset]=im_align(frames(:,:,10),ref,1);
    T2=toc;
    [T1 T2 T1/T2]
end

if 0
    %% nFrames case
    fprintf('2) nFrames case using %d frames stored on GPU already...\n',nFrames)
    % CPU
    tic
    data_matrix_cpu=zeros(nFrames,4);
    for iFrame=1:nFrames
        frame=frames(:,:,iFrame);
        [CC_max,offset]=im_align(frame,ref);
        data_matrix_cpu(iFrame,:)=[iFrame CC_max offset];
    end
    T1=toc
    
    % GPU
    tic
    g_ref=gpuArray(ref);
    g_frames=gpuArray(frames);
    data_matrix_gpu=zeros(nFrames,4,'gpuArray');
    for iFrame=1:nFrames
        g_frame=g_frames(:,:,iFrame);
        
        [CC_max,offset]=im_align(g_frame,g_ref);
        data_matrix_gpu(iFrame,:)=[iFrame CC_max offset];
    end
    T2=toc
    
    
    [T1 T2 T1/T2]
    % reset(gpuDevice(1))
    
    
    
    %%
    if all(eq(data_matrix_cpu(:),data_matrix_gpu(:)))
        disp('cool...')
    else
        %plot(data_matrix_cpu(:,2),data_matrix_gpu(:,2),'.')
    end
    
    sum(diff([data_matrix_cpu(:,2) data_matrix_gpu(:,2)],[],2).^2)
    
    
    %% convn
    A=convn(g_frame,session_data.motion_correction.kernel);
end

%% load frames onto gpu directly = slower?
nFrames=100;
tic
F1=gpuArray(session_data.get_frames(1:nFrames));
T1=toc;

tic
F2=session_data.get_frames_GPU(1:nFrames);
T2=toc;
[T1 T2 T1/T2]