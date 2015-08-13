clear all
clc

header_script

switch exp_name
    case '2015-08-10_AH02/resaved'
        iFile=22
        load_format='2015-08-10_AH02_%03d.mat';
        load_name=fullfile(data_folder,'data_analysis',sprintf('2015-08-10_AH02_%03d.mat',iFile));
        MIP_folder=fullfile(data_folder,'data_analysis','substacks');
        pixel_size_micron=[529 680]./[199 512];
        session_vector=[1 6 10 14 18 23];
    case '2015-08-10_AH03'
        iFile=1;
        load_format='2015-08-10_AH03_%03d.mat';
        load_name=fullfile(data_folder,'data_analysis',sprintf(load_format,iFile));
        MIP_folder=fullfile(data_folder,'data_analysis','substacks');
        pixel_size_micron=[500 680]./[191 512];
        session_vector=[2 6 10 14 18];
    otherwise
        die
end

%%% Revive session data for the stack session
if exist(load_name,'file')
    load(load_name,'session_data')
end

M=[cat(1,session_data.frame_info.xyz_submicron) cat(1,session_data.frame_info.laser_power)];
z_depth=M(:,3);
A=diff(z_depth)>0;
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

%% Parse out coords to fill up giant image in micron space (1 pixel is 1 micron)
coords_list=cat(1,coords.rect);
offset=min(coords_list(:,1:2));
coords_list_pos=coords_list-repmat(offset,nFiles,2);

dim=ceil(max(coords_list_pos(:,3:4)));

%% Generate blank image
STITCH.data=zeros(dim);
STITCH.gamma_val=.3;

for iFile=1:nFiles
    c=coords_list_pos(iFile,:);
    
    patch=STITCH.data(c(1)+1:c(3),c(2)+1:c(4));
    im=imresize(coords(iFile).im,size(patch));
    
    STITCH.data(c(1)+1:c(3),c(2)+1:c(4))=flipud(im);
end



%% load data sessions and plot them on top

nSessions=length(session_vector);

STITCH_add=STITCH;
for iSession=1:nSessions
    session_nr=session_vector(iSession);
    load_name=fullfile(data_folder,'data_analysis',sprintf(load_format,session_nr));
    load(load_name);
    data(iSession).center=session_data.FOV_info.center-offset([2 1]);
    %data(iSession).rect=CenterRectOnPoint([0 0 session_data.FOV_info.center([2 1])],session_data.FOV_info.center(2)-offset(2),session_data.FOV_info.center(1)-offset(1));
    data(iSession).MIP=imresize(session_data.MIP_avg.data,[382 512]);
    S(iSession)=session_data;
    
    c=data(iSession).center;
    dim=size(data(iSession).MIP);
    %[round(c(2)-dim(2)/2+1) round(c(2)+dim(2)/2) round(c(1)-dim(1)/2+1) round(c(1)+dim(1)/2)]
    STITCH_add.data(round(c(2)-dim(1)/2+1):round(c(2)+dim(1)/2),round(c(1)-dim(2)/2+1):round(c(1)+dim(2)/2))=flipud(data(iSession).MIP);
    
end
center_coords=cat(1,data.center);

%%
S.plot_FOV()

%%
session_data.imshow(STITCH,[],333)
axis xy
colormap(green)

session_data.imshow(STITCH_add)
axis xy
hold on
plot(center_coords(:,1),center_coords(:,2),'m*')
hold off

if 0
    %%
    A=S(2).MIP_avg.data;
    B=S(3).MIP_avg.data;
    session_data.imshow(A,.5,1)
    session_data.imshow(B,.5,2)
    
    [CC_max,offset]=im_align(A,B)
end