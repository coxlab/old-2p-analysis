clear all
clc

% run parse_volume_stack.m first
% works for vessel stitches as well
% make sure to add new specifications for each new session (load_format,calibration_file_name,etc.)

header_script

no_bg=2; % if 1, crop image, if 2 fit in 13x13mm space
color=1; % 1=R 2=G 3=B
offset_correction=[0 0];
save_BG_im=1;

switch exp_name
    case '2015-08-10_AH02/resaved'
        iFile_vector=22;
        load_format='2015-08-10_AH02_%03d.mat';
        calibration_file_name='AH02_20150803.mat';
        pixel_size_micron=[529 680]./[199 512];
        session_vector=[1 6 10 14 18 23];
    case '2015-08-10_AH03'
        iFile_vector=1;
        load_format='2015-08-10_AH03_%03d.mat';
        calibration_file_name='AH03_20150807.mat';
        pixel_size_micron=[500 680]./[191 512];
        session_vector=[2 6 10 14 18];
        offset_correction=[.5 -.4 0];
    case '2015-08-14_AH05'
        iFile_vector=2;
        %load_format='2015-08-14_AH05_%03d.mat';
        %calibration_file_name='AH05_20150814.mat';
        %session_vector=[3 8];
    case '2015-09-01_AH05'
        iFile_vector=1;
        pixel_size_micron=[500 680]./[191 512];
        load_format='2015-09-01_AH05_%03d.mat';
        calibration_file_name='AH05_20150901.mat';
        session_vector=[3];
        
    case '2015-08-18_AH06'
        iFile_vector=1;
        load_format='2015-08-18_AH06_%03d.mat';
        calibration_file_name='AH06_20150818.mat';
        pixel_size_micron=[500 680]./[191 512];
        session_vector=[3 6 10];
    case '2015-08-26_AH06' % red session
        iFile_vector=[1 2];
        color=2;
        load_format='2015-08-26_AH06_%03d.mat';
        calibration_file_name='AH06_20150826_v2.mat';
        pixel_size_micron=[500 680]./[191 512];
        session_vector=[];
    case '2015-09-01_AJ01' % red session
        %iFile_vector=[1 2];
        iFile_vector=3;
        color=2;
        load_format=[exp_name '_%03d.mat'];
        calibration_file_name='AJ01_20150901.mat';
        pixel_size_micron=[500 680]./[191 512];
        session_vector=[];
        
    case '082715 SR101 vessels imaging'
        iFile_vector=1:3;
        color=2;
        load_format='2015-08-27_AH04_vessels_%03d.mat';
        calibration_file_name='AH04_20150827.mat';
        pixel_size_micron=[500 680]./[191 512];
        session_vector=[];
    otherwise
        die
end

calibration_folder='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/';


%% Generate blank image
dim=[1 1]*13e3;
STITCH.data=zeros(dim);
STITCH.gamma_val=.3;

% For each file, find all coords and FOV averages
for iFile_index=1:length(iFile_vector)
    iFile=iFile_vector(iFile_index);
    load_name=fullfile(data_folder,'data_analysis',sprintf(load_format,iFile));
    %if length(iFile_vector)==1
    %    MIP_folder=fullfile(data_folder,'data_analysis','substacks');
    %else
    MIP_folder=fullfile(data_folder,'data_analysis','substacks',sprintf('session%02d',iFile));
    %end
    
    %%% Revive session data for the stack session
    if exist(load_name,'file')
        load(load_name,'session_data')
    end
    
    M=[cat(1,session_data.frame_info.xyz_submicron) cat(1,session_data.frame_info.laser_power)];
    z_depth=M(:,3);
    A=diff(z_depth)>0;
    A=medfilt1(double(A),3);
    trajectory_parts=bwlabel(A);
    trajectories=unique(trajectory_parts);
    nTrajectories=max(trajectories);
    
    session_data.get_FOV_info(pixel_size_micron)
    FOV_info=session_data.FOV_info;
    FOV_size=FOV_info.size_um;
    FOV_rect=round([0 0 FOV_size(2) FOV_size(1)]);
    
    if isdir(MIP_folder)
        files=scandir(MIP_folder,'tif');
        nFiles=length(files);
        
        coords=struct('center',[],'rect',[],'im',[]);
        for iFile=1:nFiles
            tif_name=fullfile(MIP_folder,files(iFile).name);
            im=double(imread(tif_name));
            
            %%% get position info
            idx=find(trajectory_parts==iFile);
            xy=round(M(idx(1),1:2));
            coords(iFile).center=xy;
            coords(iFile).rect=CenterRectOnPoint(FOV_rect,xy(2),xy(1));
            coords(iFile).im=im;
        end
        
    else
        disp('Folder not found')
    end
    
    %%% Parse out coords to fill up giant image in micron space (1 pixel is 1 micron)
    load(fullfile(calibration_folder,calibration_file_name))
    window_center=Calibration.window.center_coords;
    
    coords_list=cat(1,coords.rect);
    switch no_bg
        case 1
            offset=-min(coords_list(:,1:2));
            coords_list_pos=coords_list+repmat(offset,nFiles,2);
            dim=ceil(max(coords_list_pos(:,3:4)));
        case 2
            offset=round((window_center([2 1])+offset_correction([2 1]))*1000);
            dim=ceil([13 13]*1000);
            coords_list_pos=coords_list+repmat(offset,nFiles,2);
    end
    
    %%% Put FOVs in STITCH
    for iFile=1:nFiles
        c=coords_list_pos(iFile,:);
        
        patch=STITCH.data(c(1)+1:c(3),c(2)+1:c(4));
        im=imresize(coords(iFile).im,size(patch));
        
        STITCH.data(c(1)+1:c(3),c(2)+1:c(4))=flipud(im);
    end
end

%% crop STITCH
if no_bg==1
    STITCH.data=(STITCH.data(1:dim(1),1:dim(2)));
end

%%
session_data.imshow(STITCH,[],333)
hold on
circle(window_center*1000,2000,100,'r',1);
hold off
axis xy
axis equal

if color==2
    colormap(red)
end

if save_BG_im==1
    %%
    full_size=1;
    if full_size==0 
        im_folder='/Users/benvermaercke/CoxLab/MotionGUI/Images';
    else
        im_folder='/Users/benvermaercke/CoxLab/MotionGUI/Images/Full_size';
    end
    save_name=exp_name;
    save_name=fullfile(im_folder,strrep(save_name,'/','_'));
    savec(save_name)
    print(save_name,'-dpng')
    
    im=real(calc_gamma(STITCH.data,STITCH.gamma_val));
    im=im/max(im(:))*256;
    if full_size==0
        im=imresize(im,.1);
    end
    im=uint8(im);
    im=flipud(im);
    switch color
        case 1 % save data in red channel
            im=cat(3,im*0,im,im*0);
        case 2 % save data in green channel
            im=cat(3,im,im*0,im*0);
        case 3 % save data in blue channel
            im=cat(3,im*0,im*0,im);
    end
    imwrite(im,[save_name '_im.png'])
end

%% load data sessions and plot them on top
nSessions=length(session_vector);
data=struct('center',[],'MIP',[]);
if nSessions==0
else
    STITCH_add=STITCH;
    for iSession=1:nSessions
        session_nr=session_vector(iSession);
        load_name=fullfile(data_folder,'data_analysis',sprintf(load_format,session_nr));
        load(load_name);
        data(iSession).center=session_data.FOV_info.center+offset([2 1]);
        %data(iSession).rect=CenterRectOnPoint([0 0 session_data.FOV_info.center([2 1])],session_data.FOV_info.center(2)-offset(2),session_data.FOV_info.center(1)-offset(1));
        %data(iSession).MIP=imresize(session_data.MIP_avg.data,[382 512]);
        data(iSession).MIP=imresize(session_data.MIP_avg.data,[502 680]);
        %S(iSession)=session_data;
        
        c=data(iSession).center;
        dim=size(data(iSession).MIP);
        %[round(c(2)-dim(1)/2+1) round(c(2)+dim(1)/2) round(c(1)-dim(2)/2+1) round(c(1)+dim(2)/2)]
        STITCH_add.data(round(c(2)-dim(1)/2+1):round(c(2)+dim(1)/2),round(c(1)-dim(2)/2+1):round(c(1)+dim(2)/2))=flipud(data(iSession).MIP);
        
    end
    center_coords=cat(1,data.center);
    
    %%
    %S.plot_FOV()        
    session_data.imshow(STITCH_add)
    axis xy
    hold on
    plot(center_coords(:,1),center_coords(:,2),'m*')
    hold off 
    axis equal
end


if 0
    %%
    A=S(1).MIP_avg.data;
    B=S(2).MIP_avg.data;
    session_data.imshow(A,.5,1)
    session_data.imshow(B,.5,2)
    
    [CC_max,offset]=im_align(A,B)
end