clear all
clc

header_script

load_name=fullfile(data_folder,'data_analysis','2015-08-21_AH03_017.mat');

if exist(load_name,'file')==2
    load(load_name)
    
    %%
    frames=session_data.get_frames_GPU(351:400);
    
    if ismac
        %% crop frames
        frames=frames(100:150,180:240,:);
    end
    
    nRows=size(frames,1);
    nCols=size(frames,2);
    
    if ismac
        %%
        subplot(121)
        session_data.imshow(mean(frames,3))
    end
    
    %% reshape data, so we have time in rows and pixel in columns
    tic
    nFrames=size(frames,3);
    nPix=prod([nRows nCols]);
    g_M=squeeze(reshape(frames,nPix,1,nFrames))';
    size(g_M)
    toc
    
    
    %%
    if ismac
        %%
        subplot(122)
        %session_data.imshow(reshape(mean(g_M),nRows,nCols))
        session_data.imshow(g_M)
    end
    
    %% corr
    size(g_M)
    CC=corr(g_M);
    CC(between(gather(CC),[-1 1]*.4))=0;
    %A=reshape(mean(CC),nRows,nCols);
    if ismac
        imshow(abs(A),[])
    end
    
    %%
    im_save=abs(gather(A));
    im_save=im_save-min(im_save(:));
    im_save=im_save/max(im_save(:))*256;
    imwrite(uint8(im_save),'test.png')
    
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
        [ica_comps,A,W]=fastica(gather(g_M'));
        imshow(W,[])
        
    end
end
