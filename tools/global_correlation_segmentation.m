clear all
clc

header_script
load_name=fullfile(data_folder);

options.nBlocks=8;
options.frame_rate_req=4; % Hz: defines resampling rate => 4Hz = look at correlated activity in 250 ms bins
options.std_inclusion_threshold=90; % percentage
options.smoothing_px=6; % 2 for boutons, xx for cells
options.smoothing_sd=options.smoothing_px/sqrt(12); % 5 for boutons, xx for cells
options.seed_inclusion_threshold=50; % percentage
options.min_ROI_size=10;
options.save_it=0;

figure

t0=clock;
if exist(data_folder,'dir')==7
    
    expt = frGetExpt(exp_name);
    stack = readtiff(expt.dirs.regrootpn,1:options.nBlocks);
    
    %% Reduce frame_rate
    nFrames_req=32*options.frame_rate_req; % Hz
    bin_frames=round(size(stack,3)/nFrames_req)
    dec = stackGroupProject(stack,bin_frames,'sum');
    frames = single(dec);
    siz = size(dec);
    fov = imstretch(padnadapt(mean(frames,3)));
    imshow(fov)
    
    nRows=size(frames,1);
    nCols=size(frames,2);
    nFrames=size(frames,3);
    
    %% reshape data, so we have time in rows and pixel in columns
    tic
    nPix=prod([nRows nCols]);
    g_M=squeeze(reshape(frames,nPix,1,nFrames))';
    size(g_M)
    toc
    
    %% make preselection based on std per pixel
    STD=std(g_M);
    sel_activity=STD>prctile(STD,options.std_inclusion_threshold);
    tabulate(sel_activity)
    idx_list=find(sel_activity);
    
    %% corr
    g_ACT=g_M(:,sel_activity);
    CC=corr(g_ACT);
    CC_avg=mean(CC);
    CC_full=zeros(nRows*nCols,1);
    CC_full(sel_activity)=CC_avg;
    CC_2D=reshape(CC_full,nRows,nCols);
    
    
    %% smooth average CC image
    CC_avg_smooth=smooth2(CC_2D,'gauss',options.smoothing_px,options.smoothing_sd);
    
    %% take local maxima above threshold
    local_max=imregionalmax(CC_avg_smooth);
    
    if 0
        %%
        [X,Y]=find(local_max);
        imshow(calc_gamma(CC_avg_smooth,.5),[])
        hold on
        plot(Y,X,'g.')
        hold off
    end
    
    %% use those as seed points
    values=CC_avg_smooth(local_max);
    idx=find(local_max);
    [X,Y]=find(local_max);
    
    seed_points=[idx X Y values];
    sel=seed_points(:,4)>prctile(seed_points(:,4),options.seed_inclusion_threshold);
    seed_points=seed_points(sel,:);
    nSeeds=size(seed_points,1);
    
    if 0
        %%        
        imshow(calc_gamma(CC_avg_smooth,.5),[])
        hold on
        plot(seed_points(:,3),seed_points(:,2),'g.')
        hold off
    end
    
    %%        
    
    %for each seed point
    tic
    ROI_props=struct('ROI_nr',[],'idx',[],'center_x',[],'center_y',[],'CC_smoothed',[],'ROI_size',[],'CC_avg',[],'mask',[]);
    ROI_nr=1;
    for iSeed=1:nSeeds
        %% Get point index
        pxIDx=seed_points(iSeed,1);
        
        %% Grab correlation of all pixels with seed
        values=CC(idx_list==pxIDx,:);
        values(values<0.4)=0;
        if ~isempty(values)
            
            %% Reconstruct 2D image
            A=zeros(nRows,nCols);
            A(sel_activity)=values;
            %imshow(A,[])
            
            %% Select connected pixels
            SE=[0 1 0 ; 1 1 1 ; 0 1 0];
            A=imerode(A,SE);
            B=bwselect(A,seed_points(iSeed,3),seed_points(iSeed,2),8);
            A=imdilate(A,SE);
            B=imdilate(B,SE);
            ROI_size=sum(B(:));
            
            
            
            %%
            if ROI_size>options.min_ROI_size
                %% check for duplicates
                overlap=0;
                for iROI_check=1:length(ROI_props)
                    if ~isempty(ROI_props(1).ROI_nr)
                        check_01=full(ROI_props(iROI_check).mask);
                        %check_01=B;
                        check_02=B;
                        sum_im=check_01+check_02;
                        
                        overlap=2*sum(sum_im(:)==2)/(sum(check_01(:))+sum(check_02(:)));
                        %if overlap>0
                        %    imshow(sum_im,[0 2])
                        %end
                    end
                end
                
                if overlap==0
                    
                    %% Calculate stats
                    % ROI size and average correlation
                    ROI_props(ROI_nr).ROI_nr=ROI_nr;
                    ROI_props(ROI_nr).idx=seed_points(iSeed,1);
                    ROI_props(ROI_nr).center_x=seed_points(iSeed,2);
                    ROI_props(ROI_nr).center_y=seed_points(iSeed,3);
                    ROI_props(ROI_nr).CC_smoothed=seed_points(iSeed,4);
                    ROI_props(ROI_nr).ROI_size=ROI_size;
                    ROI_props(ROI_nr).CC_avg=mean(A(B==1));
                    ROI_props(ROI_nr).mask=sparse(B);
                    ROI_nr=ROI_nr+1;
                    if 0
                        %% show seed region in white and correlated regions in black
                        imshow(2*B-A,[-1 1])
                        hold on
                        plot(seed_points(iSeed,3),seed_points(iSeed,2),'ro')
                        hold off
                    end
                end
            end
        end
    end
    toc
    nROI=length(ROI_props);
    %%
    tic
    A=zeros(size(ROI_props(1).mask,1),size(ROI_props(1).mask,2),nROI);
    labeledImage=zeros(size(ROI_props(1).mask,1),size(ROI_props(1).mask,2));
    for iROI=1:nROI
        A(:,:,iROI)=full(ROI_props(iROI).mask);
        labeledImage=labeledImage+full(ROI_props(iROI).mask)*iROI;
    end
    toc
    %%
    A=mean(A,3);
    subplot(211)
    imshow(calc_gamma(CC_2D,.7),[])
    title(exp_name_txt)
    subplot(212)
    imshow(calc_gamma(A,.5),[])
    colormap(green)
    
    %% cluster based on average correlation between patches
    coloredLabels = label2rgb (labeledImage, 'cool', [1 1 1]*0 , 'shuffle');
    imshow(coloredLabels ,[])
    title(sprintf('N=%d ROIs',nROI))
    
    elapsed=etime(clock,t0);
    if options.save_it==1
        save_name=fullfile(expt.dirs.analysis,exp_name,'ROI_masks_ben.mat');
        savec(save_name)
        save(save_name,'options','nROI','ROI_props','elapsed','CC_2D','labeledImage','coloredLabels')
    end
end

