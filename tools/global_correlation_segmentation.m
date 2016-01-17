clear all
clc

header_script
load_name=fullfile(data_folder);

options.use_cache=1;
options.nBlocks=64;
options.frame_rate_req=4; % Hz: defines resampling rate => 4Hz = look at correlated activity in 250 ms bins
options.std_inclusion_threshold=90; % percentage => should become a sensible fixed value, but might depend on type of recording...
options.std_threshold=[];
options.smoothing_window_size=[11 11]; % boutons
%options.smoothing_window_size=[11 11]; % neurons
options.seed_inclusion_threshold=50; % percentage => again, we should transition to a fixed value: start recording those!
options.seed_threshold=[];
options.correlation_threshold=0.30;
options.min_ROI_size=10;
options.save_it=0;

%options.smoothing_px=4; % 2 for boutons, 6 for cells
%options.smoothing_sd=options.smoothing_px/sqrt(12); % 5 for boutons, xx for cells


figure

t0=clock;
if exist(data_folder,'dir')==7
    
    cache_name=fullfile(root_folder,'user','Ben','cache',[exp_name '.mat']);
    if options.use_cache && exist(cache_name,'file')==2
        D=load(cache_name,'options','frames');
        frames=D.frames;
        fprintf('Reloaded %d cached blocks... \n',D.options.nBlocks)
    else
        switch getenv('computername')
            case 'BKRUNCH'
                expt = frGetExpt(exp_name);
                
                filelist = dir(fullfile(expt.dirs.reggreenpn,'*.tif'));
                nFiles = length(filelist);
                
                options.nBlocks=min([options.nBlocks nFiles]);
                stack = readtiff(expt.dirs.regrootpn,1:options.nBlocks);
                
                %% Reduce frame_rate
                nFrames_req=32*options.frame_rate_req; % Hz
                bin_frames=round(size(stack,3)/nFrames_req)
                dec = stackGroupProject(stack,bin_frames,'sum');
                frames = single(dec);
            otherwise
                dataset_files=scandir(fullfile(data_folder,'data_analysis'),'mat');
                load_name=fullfile(data_folder,'data_analysis',dataset_files(3).name);
                if exist(load_name,'file')==2
                    load(load_name)
                    session_data.rebase(data_root)
                    stack=session_data.get_frames();
                end
                frames=single(stack);
        end
        %%
        if options.use_cache==1
            savec(cache_name)
            save(cache_name,'options','frames')
            disp('Done saving cache file!')
        end
    end
    
    if 0
        %%
        siz = size(dec);
        fov = imstretch(padnadapt(mean(frames,3)));
        imshow(fov)
    end
    
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
    % include temporal smoothing step to reduce noise
    STD=std(g_M);
    % include spatial smoothing step to get rid of small things getting over the
    % threshold
    options.std_threshold=prctile(STD,options.std_inclusion_threshold);
    sel_activity=STD>options.std_threshold;
    tabulate(sel_activity)
    idx_list=find(sel_activity);
    
    %% create corr matrix between all selected pixels
    g_ACT=g_M(:,sel_activity);
    CC=corr(g_ACT);
    CC_avg=mean(CC);
    CC_full=zeros(nRows*nCols,1);
    CC_full(sel_activity)=CC_avg;
    CC_2D=reshape(CC_full,nRows,nCols);
    
    
    %% smooth average CC image
    %switch getenv('computername')
    %    case 'x-BKRUNCH'
    %        CC_avg_smooth=smooth2(CC_2D,'gauss',options.smoothing_px,options.smoothing_sd);
    %    otherwise
    kernel=bellCurve2(1,options.smoothing_window_size/2,options.smoothing_window_size/6,options.smoothing_window_size,0);
    CC_avg_smooth=imfilter(CC_2D,kernel);
    %end
    
    %% take local maxima above threshold
    local_max=imregionalmax(CC_avg_smooth);
    
    if 0
        %%
        [X,Y]=find(local_max);
        
        subplot(211)
        imshow(calc_gamma(std(frames,[],3),.5),[])
        
        subplot(212)
        imshow(calc_gamma(CC_avg_smooth,.5),[])
        hold on
        plot(Y,X,'r.')
        hold off
        colormap(green)
    end
    
    %% use those as seed points
    values=CC_avg_smooth(local_max);
    idx=find(local_max);
    [X,Y]=find(local_max);
    
    seed_points=[idx X Y values];
    options.seed_threshold=prctile(seed_points(:,4),options.seed_inclusion_threshold);
    sel=seed_points(:,4)>options.seed_threshold;
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
        values(values<options.correlation_threshold)=0;
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
                    if overlap==0 && ~isempty(ROI_props(1).ROI_nr)
                        check_01=full(ROI_props(iROI_check).mask);
                        %check_01=B;
                        check_02=B;
                        sum_im=check_01+check_02;
                        
                        overlap_proportion=2*sum(sum_im(:)==2)/(sum(check_01(:))+sum(check_02(:)));
                        if overlap_proportion>0
                            overlap=1;
                            %imshow(sum_im,[0 2])
                        end
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
    nROI=length(ROI_props)
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
    
    %% Show different ROIs in different colors
    %%% X cluster based on average correlation between patches
    coloredLabels = label2rgb (labeledImage, 'cool', [1 1 1]*0 , 'shuffle');
    imshow(coloredLabels ,[])
    title(sprintf('N=%d ROIs',nROI))
    
    elapsed=etime(clock,t0);
    fprintf('Total time: %d seconds \n',round(elapsed))
    if options.save_it==1
        save_name=fullfile(expt.dirs.analysis,exp_name,'ROI_masks_ben.mat');
        savec(save_name)
        save(save_name,'options','nROI','ROI_props','elapsed','CC_2D','labeledImage','coloredLabels')
    end
    
    
    if 0
        %% Extract time traces and correlate over ROIs to merge
        trace_matrix=zeros(nFrames,nROI,'single');
        for iROI=1:nROI
            M=full(ROI_props(iROI).mask);
            sel=repmat(M,1,1,nFrames);
            R=frames.*sel;
            T=squeeze(sum(sum(R,1),2));
            trace_matrix(1:nFrames,iROI)=T;
        end
        
        %% calc correlation
        C=corr(trace_matrix);
        [C_sort,indices]=sort(median(C),'descend');
        figure(555)
        imshow(C(indices,indices),[-1 1])
        colormap jet
        
        %% cluster
        Z = linkage(C,'average');
        %[out,T]=dendrogram(Z,'ColorThreshold',3)
        cluster_allocation_vector=cluster(Z,'CutOff',1+3/100)
        %cluster_allocation_vector=cluster(Z,'MaxClust',27,'Criterion','distance')
        nClusters=length(unique(cluster_allocation_vector))
        
        %%% show spatially clustered ROIs
        clusterImage=zeros(size(ROI_props(1).mask,1),size(ROI_props(1).mask,2));
        for iROI=1:nROI            
            cluster_nr=cluster_allocation_vector(iROI);
            clusterImage=clusterImage+full(ROI_props(iROI).mask)*cluster_nr;
        end
        coloredLabels = label2rgb (clusterImage, 'cool', [1 1 1]*0 , 'shuffle');
        imshow(coloredLabels ,[])
        title(sprintf('N=%d clusters',nClusters))
        
    end
end

