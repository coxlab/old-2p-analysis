clear all
clc

header_script


if ismac
    use_GPU=0;
else
    use_GPU=0;
end
load_name=fullfile(data_folder);

if exist(data_folder,'dir')==7
    
    expt = frGetExpt(exp_name);
    stack = readtiff(expt.dirs.regrootpn,1:65);
    %dec = stackGroupProject(stack,256,'sum');
    
    %%
    frame_rate_req=4; % Hz
    nFrames_req=32*frame_rate_req; % Hz
    bin_frames=round(size(stack,3)/nFrames_req);
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
    
    sel_activity=STD>prctile(STD,96);
    tabulate(sel_activity)
    %sel_activity=STD>20;
    if ismac
        %%
        %hist(gather(STD),length(STD)/20)
        figure(444)
        subplot(121)
        session_data.imshow(reshape(sel_activity,nRows,nCols),.7)
    end
    
    %%
    if ismac
        %%
        %subplot(122)
        %session_data.imshow(reshape(mean(g_M),nRows,nCols))
        %session_data.imshow(g_M)
    else
        
    end
    
    %% corr
    g_ACT=g_M(:,sel_activity);
    size(g_ACT)
    CC=corr(g_ACT);
    %CC(between(gather(CC),[-1 1]*.3))=0;
    CC_avg=mean(CC);
    if use_GPU==1
        CC_full=zeros(nRows*nCols,1,'gpuArray');
    else
        CC_full=zeros(nRows*nCols,1);
    end
    CC_full(sel_activity)=CC_avg;
    A=reshape(CC_full,nRows,nCols);
    
    g_ACT_CORR=g_ACT(:,CC_avg>.01);
    size(g_ACT_CORR)
    if ismac
        %imshow(abs(A),[])
    end
    
    %%
    if use_GPU==1
        im_save=abs(gather(A));
    else
        im_save=abs(A);
    end
    im_save=im_save-min(im_save(:));
    im_save=im_save/max(im_save(:))*256;
    save_name=fullfile(data_folder,'test.png');
    %imwrite(uint8(im_save),save_name)
    
    if ismac
        %%
        %subplot(122)
        imshow(im_save,[0 256])               
        
        %click a pixel and get all pixels that correlate with this one
        [c,r,button]=ginput(1);
        
        %% convert to pix index
        %pxIDx=round((c-1)*nRows+r);
        pxIDx=find(im_save(:)==max(im_save(:)),1,'first'); % choose max corr pix
        
        idx_list=find(sel_activity);
        values=CC(idx_list==pxIDx,:);
        
        A=im_save;
        TH=7/10;
        values(values<TH)=0;
        A(sel_activity)=values;
        
        imshow(A,[-1 1])
        hold on
        plot(c,r,'r*')
        hold off
        axis ij
        %session_data.imshow(im_save)
        
        %% smooth average CC image
        CC_avg_smooth=smooth2(im_save,'gauss',5,2);
        imshow(calc_gamma(CC_avg_smooth,.5),[])
        
        %% take local maxima above threshold
        local_max=imregionalmax(CC_avg_smooth);
        imshow(local_max)
        
        %% use those as seed points
        values=CC_avg_smooth(local_max);
        idx=find(local_max);
        [X,Y]=find(local_max);
        
        seed_points=[idx X Y values];
        sel=seed_points(:,4)>prctile(seed_points(:,4),80);
        seed_points=seed_points(sel,:);
        nSeeds=size(seed_points,1);
        
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
            if sum(values)>0
                
                %% Reconstruct 2D image
                A=zeros(nRows,nCols);
                A(sel_activity)=values;
                
                %% Select connected pixels
                SE=[0 1 0 ; 1 1 1 ; 0 1 0];
                A=imerode(A,SE);
                B=bwselect(A,seed_points(iSeed,3),seed_points(iSeed,2),8);
                A=imdilate(A,SE);
                B=imdilate(B,SE);
                ROI_size=sum(B(:));
                
                if ROI_size>5
                    %% Calculate stats
                    % ROI size and average correlation
                    %[sum(B(:)) mean(A(B==1))]
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
        toc
        
        nROI=length(ROI_props);
        
        %%  
        tic
        A=zeros(size(ROI_props(1).mask,1),size(ROI_props(1).mask,2),nROI);
        for iROI=1:nROI
            A(:,:,iROI)=full(ROI_props(iROI).mask);            
        end
        toc
        A=mean(A,3);
        imshow(A,[])
        
        %% cluster based on average correlation between patches
        
    end
    die
    
    %% create distance matrix
    D=1-CC;
    % 1+1-R
    
    
    mappedX=mds(D,5);
    
    
    nClusters=5;
    cluster_mapping=kmeans(gather(mappedX),nClusters);
    tabulate(cluster_mapping)
    
    if ismac
        colors={'r','g','b'};
        clf
        hold on
        for iClust=1:nClusters
            sel=cluster_mapping==iClust;
            plot(mappedX(sel,1),mappedX(sel,2),'*','color',colors{iClust})
        end
        hold off
        axis equal
        
        %% show clustered image
        clf
        im=reshape(cluster_mapping,nRows,nCols);
        imagesc(im)
        die
    end
    
    
    %%
    
    if ismac&&0
        %%
        session_data.imshow(CC,1.2)
        axis square
        
        %%
        %Z=linkage(center_coords,'average');
        %C=cluster(Z,'cutoff',min_dist,'Criterion','distance');
        [coeff, score, latent, tsquared, explained, mu] = pca(g_M);
        
    end
    
    if 0
        
        
        %% linkage
        %Z=linkage(gather(g_M'));
        %dendrogram(Z)
        
        %% ICA
        [ica_comps,A,W]=fastica(gather(g_ACT_CORR'));
        imshow(W,[])
        
    end
end
