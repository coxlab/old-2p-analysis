clear all
clc

%%% point to FOV-file
FOV_file_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-10_AF11_exp/data_analysis/20150410_AF11_004.mat';

%%% point to stack-file
stack_file_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-10_AF11_exp/20150410_AF11_005.tif';


%%% Load MIP for the selected FOV
load(FOV_file_name,'session_data')
FOV=session_data.MIP_avg.data;
%%

%%% Load in stack
%file_name=fullfile(data_folder,files(iFile).name);
info=imfinfo(stack_file_name);

%%% Get basic info about movie
nFrames=length(info);
W=info(1).Width;
H=info(1).Height;

%%% Read stack
offset=double(intmax('uint16')/2);
frame_data=zeros(nFrames,10);
frames=zeros(H-1,W,nFrames);
for iFrame=1:nFrames
    data=double(imread(stack_file_name,iFrame,'info',info)); % no flyback line
    frame=data(1:end-1,:);
    flyBack=data(end,end-17:end);
    
    date_num=flyBack(1:6)./[1 1e3 1e3 1e3 1e3 1e3];
    xyz_hi_res=(flyBack(11:13)-offset)/10; % micron
    xyz=flyBack(14:16)-offset; % micron
    piezo=flyBack(17); % micron
    laser_power=flyBack(18); % percent! not mW
    
    frame_data(iFrame,:)=[iFrame datenum(date_num) xyz xyz_hi_res piezo laser_power];
    frames(:,:,iFrame)=frame;
end

%%
dataMatrix=zeros(nFrames,5);
t0=clock;
for iFrame=1:nFrames
    frame=frames(:,:,iFrame);
    [CC_max,offset]=im_align(frame,FOV);
    dataMatrix(iFrame,:)=[iFrame mean(frame(:)) CC_max offset];
    progress(iFrame,nFrames,t0)
end

%% calc show frames
frames_show=frames;
frames_show=calc_gamma(frames_show,gamma_val);
frames_show=medfilt3(frames_show,[1 1 3]);


%% show data
TH=.60; % minimum correlation to reject matches
Z_offset=0; % where does surface start?
gamma_val=.5;
green=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];

CC=dataMatrix(:,3);


FOV_depth=mode(session_data.dataMatrix(:,5));


%true_depth=FOV_depth-Z_axis(1);
true_depth=find(frame_data(:,8)==FOV_depth)
time_axis=frame_data(:,1);
Z_axis=frame_data(:,8);
Z_axis=Z_axis-Z_axis(1)+Z_offset;
laser_power=frame_data(:,10);
[m,loc]=max(CC);


clf
if m<TH
    error('No match was found...')
else
    switch 2
        case 1
            match_frame=frames(:,:,loc);
        case 2
            match_frame=mean(frames(:,:,loc-1:loc+1),3);
        case 3
            match_frame=medfilt3(frames(:,:,loc-1:loc+1),[1 1 1]);
    end
    
    ROI=session_data.ROI_definitions(2).ROI;
    nROI=length(ROI);
    
    figure(1)
    
    
    %subplot(223)
    %imshow(calc_gamma(FOV,gamma_val),[])
    
    subplot(6,1,[1 5])
    h=imshow(calc_gamma(match_frame,gamma_val),[]);
    hold on
    for iROI=1:nROI
        plot(ROI(iROI).coords_MIP(:,1)+dataMatrix(loc,4),ROI(iROI).coords_MIP(:,2)+dataMatrix(loc,5),'color',[1 1 1]*.7)
    end
    hold off
    title(sprintf('%d and %d',dataMatrix(loc,4:5)))
    colormap(green)
        
    subplot(616)
    plot(time_axis,CC)
    hold on
    plot(time_axis([loc loc]),[0 max(CC)*1.2],'r-')
    plot(time_axis,laser_power/max(laser_power),'r')
    plot(time_axis,Z_axis/max(Z_axis),'g')
    plot([true_depth true_depth],[0 max(CC)*1.2],'c')
    %text(Z_axis(loc),max(CC)*.2,sprintf('%3.1f',Z_axis(loc)))
    indicator=plot(time_axis([loc loc]),[0 max(CC)*1.2],'k-');
    hold off
    box off
    axis([0 nFrames 0 1])
    
    set(gcf,'WindowButtonMotionFcn',{@moveFcn_stack,indicator,h,frames_show,frame_data})
end

