clear all
clc

%%% Create a more elegant, automated way to create the illustrator files
%%% Find a way to get information about surgery image, injection sites and
%%% epi images to come together, the same frame should hold 2p images we
%%% may have stitched together and guide further action.

%%% This would entail an animal specific file or a file/database with
%%% records for all animals.
% Fields:
% - surgery image (+scaling factor) : arbitrary units, depend on zoomfactor
% - injection coordinates (+offset) => mm stereotax
% - epi images (+scaling factor and offset) => mm xyz-stage
% - blood vessel SR101 (+scaling factor and offset)  => mm xyz-stage
% - 2p images  (+scaling factor and offset)  => mm xyz-stage
% - window coordinates (3D)  => mm xyz-stage

% This could form the basis of various stitching paths

% We need MotionGUI to drive our camera based acquisition of epi stitches.
% Based on FOV size in microns, move a certain percentage in all directions
% and capture new image, naming can be manual.

animal_selector=4;
animal_list={'AF11','AF17','AG01','AG02'};
animal_ID=animal_list{animal_selector};
folder_name='2015-06-18_AF17_init';

im_root_folder=sprintf('/Users/benvermaercke/Dropbox (coxlab)/2p-data/surgery_rig_images/%s',animal_ID);
epi_root_folder=sprintf('/Users/benvermaercke/Dropbox (coxlab)/2p-data/%s',folder_name);


switch animal_ID
    case 'AF11'
        craniotomy_center=[-7.5 4]; % AP ML in mm, recorded during implant. Final center may be off by exact placement of headplate and location of drilling.
        
        surgery_im_name='IMG_5675.JPG'; % will be converted into a micron sized grid
        real_diameter=4; % mm, approximate coverslip diameter
        surgery_im_scaling_factor=343.5383; % how many pixels for 1 mm
        surgery_im_center=[1553 1273]; % in pixels
        orientation=1;
        
        injection_coords_offset=[-1.2 -.8]; % in mm, how we have to shift to coregister injections with surgery image
        injection_coords=[0 0 ; .43 .65 ; .48 1.21 ; 1.14 1.72 ; .45 1.97]; % AP ML in mm, recorded during window surgery
        injection_speed=50; % nl/min
        injection_volumes=[500 250 250 250 250]; % nl
        injection_depths=[-.62 -.65 -.64 -.62 -.72]; % mm, from pial surface
        injection_qualities=[3 1 3 4 5]; % scale of 1 to 5
        injection_names={'IMG_5342.jpg','','IMG_5345.jpg','IMG_5350.jpg','IMG_5354.jpg'};
        
        epi_names={'01_V1_medial.jpg','02_V1_central.jpg','03_V1_medial_lateral.jpg','04_V1_lateral.jpg','05_V1_anterior_lateral.jpg'};
        epi_FOV_size_micron=[945.9854 709.4891]/1000; % in mm
        epi_FOV_size_px=[640 486]; % px
        epi_scaling_factor=epi_FOV_size_px/epi_FOV_size_micron;
        epi_coordinates=[ -1.3990 -0.8651 0.7748 ; -0.5783 -0.2790 0.7748 ; 0.0232 -0.0488 0.7034 ; 0.8637 -0.2745 0.7034 ; 0.5270 0.4743 0.7034]; % in mm
        epi_offset=[0 0]; % mm
        epi_coordinates_px=epi_coordinates(:,1:2)*surgery_im_scaling_factor+repmat(surgery_im_center+epi_offset*surgery_im_scaling_factor,size(epi_coordinates,1),1);
        
    case 'AF17'
        craniotomy_center=[-8.49 4.5]; % AP ML in mm, recorded during implant. Final center may be off by exact placement of headplate and location of drilling.
        
        real_diameter=4; % mm, approximate coverslip diameter
        switch 2
            case 1
                surgery_im_name='IMG_6275.jpg'; % will be converted into a micron sized grid
                orientation=2;
                surgery_im_scaling_factor=418.322; % how many pixels for 1 mm
                surgery_im_center=[1598 1291]; % in pixels
                injection_coords_offset=[-.95 -.8]; % in mm, how we have to shift to coregister injections with surgery image
            case 2
                surgery_im_name='IMG_6479.JPG'; % will be converted into a micron sized grid
                orientation=1;
                surgery_im_scaling_factor=355.347;
                surgery_im_center=[1791 1392];
                injection_coords_offset=[-.95 -.8]; % in mm, how we have to shift to coregister injections with surgery image
        end
        
        
        
        injection_coords=[0 0 ; .02 .40 ; .42 .48 ; -.29 .91 ; -.25 1.56 ; -.12 2.42 ; 1.36 2.51 ; .74 .71 ]; % AP ML in mm, recorded during window surgery
        injection_speed=50; % nl/min
        injection_volumes=[150 150 150 150 150 150 150 150]; % nl
        injection_depths=[-.58 -.57 -.62 -.63 -.61 -.60 -.59 -.60]; % mm, from pial surface
        injection_qualities=[3 3 1 1 4 4 5 3]; % scale of 1 to 5
        injection_names={'IMG_6086.jpg','IMG_6090.jpg','IMG_6092.jpg','IMG_6094.jpg','IMG_6102.jpg','IMG_6107.jpg','IMG_6109.jpg',''};
        
        epi_names={'center_BF.jpg','center_posterior_BF.jpg'};
        epi_FOV_size_micron=[945.9854 709.4891]/1000; % in mm
        epi_FOV_size_px=[640 486]; % px
        epi_scaling_factor=epi_FOV_size_px/epi_FOV_size_micron;
        epi_coordinates=[0 0 .56 ; -0.2621 -0.6227 0.8216]; % in mm
        epi_offset=[.2 -.2]; % mm
        epi_coordinates_px=epi_coordinates(:,1:2)*surgery_im_scaling_factor+repmat(surgery_im_center+epi_offset*surgery_im_scaling_factor,size(epi_coordinates,1),1);
        
    case 'AG01'
        craniotomy_center=[-7.89 5]; % AP ML in mm, recorded during implant. Final center may be off by exact placement of headplate and location of drilling.
        
        surgery_im_name='IMG_6491.JPG'; % will be converted into a micron sized grid
        real_diameter=4; % mm, approximate coverslip diameter
        
        surgery_im_scaling_factor=348.376; % how many pixels for 1 mm
        surgery_im_center=[1546 1377]; % in pixels
        
        orientation=1;
        
        injection_coords_offset=[-.9 -.1]; % in mm, how we have to shift to coregister injections with surgery image
        injection_coords=[0 0 ; -.71 -.01 ; -.57 .87 ; -.87 1.20 ; -.84 1.78 ; -.44 2.08 ; .12 2.31 ; .81 1.36 ;]; % AP ML in mm, recorded during window surgery
        injection_speed=50; % nl/min
        injection_volumes=[150 150 150 150 150 150 150 150]; % nl
        
        injection_depths=[-.60 -.60 -.60 -.60 -.60 -.60 -.60 -.60]; % mm, from pial surface
        injection_qualities=[1 1 3 1 1 1 1 1]; % scale of 1 to 5
        injection_names={'IMG_6180.jpg','IMG_6184.jpg','IMG_6187.jpg','IMG_6192.jpg','IMG_6196.jpg','IMG_6200.jpg','IMG_6205.jpg','IMG_6209.jpg'};
        
        epi_names={''};
        epi_FOV_size_micron=[945.9854 709.4891]/1000; % in mm
        epi_FOV_size_px=[640 486]; % px
        epi_scaling_factor=epi_FOV_size_px/epi_FOV_size_micron;
        epi_coordinates=[ 0 0 ]; % in mm
        epi_offset=[0 0]; % mm
        epi_coordinates_px=epi_coordinates(:,1:2)*surgery_im_scaling_factor+repmat(surgery_im_center+epi_offset*surgery_im_scaling_factor,size(epi_coordinates,1),1);
    case 'AG02'
        craniotomy_center=[-7.99 5]; % AP ML in mm, recorded during implant. Final center may be off by exact placement of headplate and location of drilling.
        
        surgery_im_name='IMG_6497.JPG'; % will be converted into a micron sized grid
        real_diameter=4; % mm, approximate coverslip diameter
        
        surgery_im_scaling_factor=378.649; % how many pixels for 1 mm
        surgery_im_center=[1455 1334]; % in pixels
        
        orientation=1;
        
        injection_coords_offset=[-1.5 -0]; % in mm, how we have to shift to coregister injections with surgery image
        injection_coords=[0 0 ; .43 .82 ; 1.16 1.18 ; -1.21 1.43 ; -1.04 2.26 ; -1.02 2.74 ; -.15 2.56 ; -.52 3.17 ; .41 2.39]; % AP ML in mm, recorded during window surgery
        injection_speed=50; % nl/min
        injection_volumes=[150 150 150 150 150 150 150 150 150]; % nl
        
        injection_depths=[-.60 -.60 -.85 -.60 -.60 -.60 -.60 -.60 -.60]; % mm, from pial surface
        injection_qualities=[1 1 1 1 3 1 2 1 3]; % scale of 1 to 5
        injection_names={'IMG_6226.jpg','IMG_6229.jpg','IMG_6231.jpg','IMG_6235.jpg','IMG_6237.jpg','IMG_6245.jpg','IMG_6255.jpg','IMG_6259.JPG',''};
        
        epi_names={''};
        epi_FOV_size_micron=[945.9854 709.4891]/1000; % in mm
        epi_FOV_size_px=[640 486]; % px
        epi_scaling_factor=epi_FOV_size_px/epi_FOV_size_micron;
        epi_coordinates=[ 0 0 ]; % in mm
        epi_offset=[0 0]; % mm
        epi_coordinates_px=epi_coordinates(:,1:2)*surgery_im_scaling_factor+repmat(surgery_im_center+epi_offset*surgery_im_scaling_factor,size(epi_coordinates,1),1);
        
    otherwise
        disp('No data for this rat')
end

% resize image so 1 pixel is 1 micron
if 0 % if you need help quantifying these dimensions
    %%
    loadName=fullfile(im_root_folder,surgery_im_name);
    if orientation==1
        BG=fliplr(imread(loadName));
    else
        BG=flipud(imread(loadName));
    end
    
    warning off
    imshow(BG,[])
    axis xy
    axis xy
    warning off
    hold on
    p=plot(0,0,'m*');
    hold off
    title('mark 4 or more point on the edge of the coverslip')
    
    %%% Collect calibration points
    coords=[];
    running=1;
    while running==1
        [x,y,button]=ginput(1);
        if button==1
            coords=cat(1,coords,[x y]);
            set(p,'xData',coords(:,1),'yData',coords(:,2))
        else
            running=0;
        end
    end
    
    ellipse=fit_ellipse(coords(:,1),coords(:,2));
    center=[ellipse.X0_in ellipse.Y0_in 0];
    diameter=mean([ellipse.short_axis ellipse.long_axis]);
    radius=diameter/2;
    hold on
    plot(center(1),center(2),'ro')
    plotCircle(center,radius,100,'r-');
    hold off
    axis equal
    
    surgery_im_scaling_factor=diameter/real_diameter;
    
    sprintf('surgery_im_scaling_factor=%3.3f;\nsurgery_im_center=[%d %d];\n',[surgery_im_scaling_factor round(center(1:2))])
end


% concat injections
nInjections=size(injection_coords,1);
injections=struct;
for iInject=1:nInjections
    injections(iInject).coords_mm=fliplr(injection_coords(iInject,:));
    injections(iInject).coords_px=injections(iInject).coords_mm*surgery_im_scaling_factor+surgery_im_center+injection_coords_offset*surgery_im_scaling_factor;
    injections(iInject).speed=injection_speed;
    injections(iInject).volume=injection_volumes(iInject);
    injections(iInject).depth=injection_depths(iInject);
    injections(iInject).quality=injection_qualities(iInject);
    injections(iInject).image_name=injection_names{iInject};
end
injection_sites=cat(1,injections.coords_px);
injection_sites_mm=cat(1,injections.coords_mm);
injection_spacing=calc_dist([injection_sites_mm(1:end-1,:) injection_sites_mm(2:end,:)])*1e3;


% concat epi fields
nEpi=length(epi_names);
epi=struct;
for iEpi=1:nEpi
    epi(iEpi).im_name=epi_names{iEpi};
    epi(iEpi).coords=epi_coordinates(iEpi,:);
    %epi(iEpi).coords_px=epi(iEpi).coords(1:2)*surgery_im_scaling_factor+surgery_im_center;
    epi(iEpi).coords_px=epi_coordinates_px(iEpi,1:2);
    
    %epi(iEpi).im_name=epi_names{iEpi};
    %epi(iEpi).im_name=epi_names{iEpi};
end
epi_locations=cat(1,epi.coords_px);


%% %% plot final product

loadName=fullfile(im_root_folder,surgery_im_name);
if orientation==1
    BG=fliplr(imread(loadName));
else
    BG=flipud(imread(loadName));
end

figure(animal_selector)
imshow(BG,[])
axis xy

hold on
plotCircle([surgery_im_center 0],real_diameter/2*surgery_im_scaling_factor,100,'r-');
plot(surgery_im_center(1),surgery_im_center(2),'w+')
plot(injection_sites(:,1),injection_sites(:,2),'rv')

for iInj=1:nInjections-1
    plot([injection_sites(iInj,1) injection_sites(iInj+1,1)],[injection_sites(iInj,2) injection_sites(iInj+1,2)],'r')
    t=text(mean([injection_sites(iInj,1) injection_sites(iInj+1,1)]),mean([injection_sites(iInj,2) injection_sites(iInj+1,2)]),sprintf('%3d',round(injection_spacing(iInj))));
    set(t,'color','w')
end
%plot(epi_coordinates_px(:,1),epi_coordinates_px(:,2),'ws')

%plot(epi_locations(:,1),epi_locations(:,2),'ks')
hold off
axis equal
title(animal_ID)


%%% make a brightfield stitch as a base, dimensions are off using
%%% dissection scope