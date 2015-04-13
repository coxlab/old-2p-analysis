clear all
clc

[~, user_name] = system('whoami');user_name=char(user_name);
root_folder=fullfile('/Users',user_name,'/Dropbox (coxlab)');
% define data folder
data_root=fullfile(root_folder,'2p-data');

data_folder=fullfile(data_root,'2015-03-03_AF03-light_awake');
file_name='20150303_AF03_001.tif';

loadName=fullfile(data_folder,file_name);
if exist(loadName,'file')==2
    
    tic
    info=imfinfo(loadName);
    nFrames=length(info);
    W=info.Width;
    H=info.Height;
    
    %nFrames=min([1000 nFrames]);
    
    frames=zeros(H-1,W,nFrames);
    for iFrame=1:nFrames
        frame=double(imread(loadName,iFrame,'info',info));
        frames(:,:,iFrame)=double(frame(1:end-1,:)); % no flyback line double
    end
    fprintf('Reading file took: %3.2fs\n',toc)
end

%% Get average or std of N frames
tic
block_size=20;
nBlocks=nFrames/block_size;
blocks=zeros(H-1,W,nBlocks);
for iBlock=1:nBlocks
    sel=(iBlock-1)*block_size+1:iBlock*block_size;
    blocks(:,:,iBlock)=std(frames(:,:,sel),[],3);
end
fprintf('Calculating block averages took: %3.2fs\n',toc)

%% Calculate CC between succesive frames
tic
T=(1:nBlocks-1)*block_size;
data_matrix=zeros(nBlocks-1,8);
for iBlock=1:nBlocks-1
    im1=blocks(:,:,iBlock);
    im2=blocks(:,:,iBlock+1);
    CC=normxcorr2(im2,im1);
    [x,y]=find(CC==max(CC(:)));
    x_shift=x-size(im1,1)/2+1;
    y_shift=y-size(im1,2)/2+1;
    data_matrix(iBlock,:)=[iBlock iBlock+1 T(iBlock) max(CC(:)) x y x_shift y_shift];
end
fprintf('Calculating offsets took: %3.2fs\n',toc)

%%

plot(data_matrix(:,8),data_matrix(:,7),'*')

axis([0 H 0 W])