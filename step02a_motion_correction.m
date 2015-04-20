%%% BV20150410: add in an coarser timeline, averaging n frames and check
%%% the shift over those. This should be more robust and could serve as an
%%% approximate ground truth for the frame by frame analysis
% just correcting x-y translation for now. angles do not seem to very at
% the timescale of 1 session

%%% BV20150410: consider a hybrid approach where we use the coarse
%%% alignment first and then get offset from each single frame to the
%%% coarse image.
%%% Two advantage, stop large drifts based on minor shifts due to error and
%%% more rigorous comparison because 1 of the two images in the fine
%%% comparison is more reliable. Problem arises when large shift renders
%%% the average image useless.
%%% For now, doing just the fine alignment is not entirely robust

clear all
clc

apply_motion_correction=0;
save_it=1;

header_script

%%
%%% Use uigetdir for general use
cd(data_root)
data_folder=uigetdir(data_root);

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
    
    tic
    frames=zeros(rows-1,cols,nFrames);
    for iFrame=1:nFrames
        frame=double(imread(session_data.file_name,iFrame,'info',info));
        frames(:,:,iFrame)=double(frame(1:end-1,:)); % no flyback line double
    end
    fprintf('Loading frames took %3.2f seconds.\n',toc)
    
     %% get number of blank frame to ignore at start
     mean_lum=squeeze(mean(mean(frames,1),2));
     blank_frames=false(size(mean_lum));
     blank_frames(1:find(zscore(mean_lum(1:round(end/4)))<-4,1,'last')+1)=true;
     
    
    if apply_motion_correction==1
        %% create block averages over n frames
        tic
        block_size=15;
        nBlocks=floor(nFrames/block_size);
        blocks=zeros(rows-1,cols,nBlocks);
        for iBlock=1:nBlocks
            sel=(iBlock-1)*block_size+1:iBlock*block_size;
            blocks(:,:,iBlock)=mean(frames(:,:,sel),3);
        end
        fprintf('Calculating block averages took: %3.2fs\n',toc)
        
        %%
        gamma_val=.5;
        shift_matrix=zeros(nFrames,5);
        motion_correction.kernel_size=[6 6];
        motion_correction.sigma=.65;
        kernel=bellCurve2(1,motion_correction.kernel_size/2,[motion_correction.sigma motion_correction.sigma],motion_correction.kernel_size,0);
        motion_correction.kernel_size_coarse=[6 6];
        motion_correction.sigma_coarse=1.1;
        kernel_coarse=bellCurve2(1,motion_correction.kernel_size_coarse/2,[motion_correction.sigma_coarse motion_correction.sigma_coarse],motion_correction.kernel_size_coarse,0);
        
        motion_correction.downsample_factor=2;
        motion_correction.correction_method=2;
        
        %% do coarse analysis
        tic
        T=(1:nBlocks)*block_size;
        data_matrix=zeros(nBlocks-1,8);
        data_matrix(1,:)=[0 1 T(1) 0 0 0 0 0];
        for iBlock=2:nBlocks
            im1=blocks(:,:,iBlock-1);
            im2=blocks(:,:,iBlock);
            
            %%% Gaussian filter image with .65px kernel
            im1_smooth=convn(im1,kernel_coarse,'same');
            im2_smooth=convn(im2,kernel_coarse,'same');
            
            %%% Resize image 2x or 4x
            resize_factor=1/motion_correction.downsample_factor;
            im1_smooth_resize=imresize(im1_smooth,resize_factor);
            im2_smooth_resize=imresize(im2_smooth,resize_factor);
            
            CC=normxcorr2(im1_smooth_resize,im2_smooth_resize);
            [x,y]=find(CC==max(CC(:)));
            x_shift=x-size(im1_smooth_resize,1);
            y_shift=y-size(im1_smooth_resize,2);
            fprintf('%d %d\n',[x_shift y_shift])
            data_matrix(iBlock,:)=[iBlock-1 iBlock T(iBlock) max(CC(:)) x_shift/resize_factor y_shift/resize_factor 0 0];
        end
        data_matrix(:,7)=cumsum(data_matrix(:,5));
        data_matrix(:,8)=cumsum(data_matrix(:,6));
        data_matrix
        fprintf('Calculating offsets took: %3.2fs\n',toc)
        
        %% do fine analysis
        t0=clock;
        tic
        shift_matrix(1,:)=[1 0 0 0 0];
        for iFrame=2:nFrames
            if ismember(iFrame,find(blank_frames))
                shift_matrix(iFrame,:)=[iFrame 0 0 0 0];
            else
                % Grab frame for processing
                F1=frames(:,:,iFrame-1);
                F2=frames(:,:,iFrame);
                
                %%% Gaussian filter image with .65px kernel
                F1_smooth=convn(F1,kernel,'same');
                F2_smooth=convn(F2,kernel,'same');
                
                %%% Resize image 2x or 4x
                resize_factor=1/motion_correction.downsample_factor;
                F1_smooth_resize=imresize(F1_smooth,resize_factor);
                F2_smooth_resize=imresize(F2_smooth,resize_factor);
                
                %%% Estimate displacement: has to be spot on because error will
                %%% accumulate and can cause a drift.
                A=F1_smooth_resize;
                B=F2_smooth_resize;
                
                switch motion_correction.correction_method
                    case 1 % normxcorr2, more accurate but slower
                        cc=normxcorr2(A,B);
                        [peakx,peaky]=find(cc==max(cc(:)));
                        r=peakx-size(A,1);
                        c=peaky-size(A,2);
                        
                    case 2 % fft phase matching: fast and often good enough
                        [r,c]=PCdemo(A,B);
                end
                
                shift_matrix(iFrame,:)=[iFrame r/resize_factor c/resize_factor 0 0];
                progress(iFrame,nFrames,t0)
                
                if 0
                    %%
                    subplot(121)
                    imshow(calc_gamma(F1_smooth,gamma_val),[])
                    hold on
                    plot(size(F1_smooth,1)/2,size(F1_smooth,2)/2,'m*')
                    hold off
                    subplot(122)
                    imshow(calc_gamma(F2_smooth,gamma_val),[])
                    hold on
                    plot(size(F2_smooth,1)/2+r,size(F2_smooth,2)/2+c,'m*')
                    hold off
                    title([r c])
                    colormap(green)
                end
            end
        end
        shift_matrix(:,4)=cumsum(shift_matrix(:,2));
        shift_matrix(:,5)=cumsum(shift_matrix(:,3));
        toc
        
        
        %% check quality of the correction
        for iFrame=1:nFrames
            iCoarse=find(T<iFrame,1,'last');
            if isempty(iCoarse)
                %D(iFrame)=1;
                A=data_matrix(1,7:8);
            else
                %D(iFrame)=iCoarse;
                A=data_matrix(iCoarse,7:8);
            end
            B=shift_matrix(iFrame,4:5);
            d(iFrame)=calc_dist([A B]);
        end
        if sum(d)>100
            plot(d)
            error('Too much difference between coarse and fine motion corrections...')
            
            shift_matrix(:,2:5)=data_matrix(5:8)
        end
        
        %% visualize
        plot(shift_matrix(:,4))
        hold on
        plot(shift_matrix(:,5),'r')
        plot(data_matrix(:,3),data_matrix(:,7),'g')
        plot(data_matrix(:,3),data_matrix(:,8),'c')
        hold off
        motion_correction.shift_matrix=shift_matrix;
        
        if 0
            %% show movie
            for iFrame=1:nFrames
                iCoarse=find(T<iFrame,1,'last');
                if isempty(iCoarse)
                    D(iFrame)=1;
                else
                    D(iFrame)=iCoarse;
                end
            end
            
            %%
            for iFrame=1:nFrames
                offset=[size(frames,1)/2 size(frames,2)/2]
                imshow(calc_gamma(frames(:,:,iFrame),gamma_val),[])
                hold on
                plot(offset(1)+shift_matrix(iFrame,5),offset(2)+shift_matrix(iFrame,4),'m*')
                plot(offset(1)+data_matrix(D(iFrame),8),offset(2)+data_matrix(D(iFrame),7),'mo')
                hold off
                colormap(green)
                drawnow
            end
        end
    end
    %%
    
    if save_it==1
        %%
        session_data.mean_lum=mean_lum;
        session_data.blank_frames=blank_frames;
        if apply_motion_correction==1
            session_data.motion_correction=motion_correction;
        end
        save(loadName,'session_data');
    end
end