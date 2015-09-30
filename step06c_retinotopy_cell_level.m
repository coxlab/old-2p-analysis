% Cool script so far, outlines where we need to go but is still quite slow
% and inefficient. We need to get into a scenario where we have a flat
% structure of all cells, where we can add data as data comes in. So all
% data can be plotted in one go. Analysis happens one time when ROIs are
% defined. Every ROI has links to Animal, Session and FOV with appropriate
% calibration file/coordinate system used on for that session. So each ROI
% can be placed in absolute space without issues. 

% ROI has absolute coordinates (I smell a class coming), 40x40 mask around
% that coordinate, reference to animal, session, FOV, exp_type...
% each session will yield a file that contains all ROIs of type ROI_class
% these can be stored in the same folder and data can be requested per
% animal, all session or selected, all FOVs or selected.

% Need a vessel image to more accurately verify FOV position within window.

clear all
clc

header_script
animal_ID='AH03';

offset_correction=[0 0];
switch animal_ID
    case 'AH02'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH02_20150803.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH02_resaved_im.png';
    case 'AH03'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH03_20150807.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH03_im.png';
        offset_correction=[-.4 .5];
    case 'AH05'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH05_20150901.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-09-01_AH05_im.png';
end

dataset_folder=fullfile(dataset_root,animal_ID);

%%% First try this: on the background image, plot all the ROIs
% this will require to have the center of each FOV and then plot each ROI,
% defined relative to the FOV in the coordinate system of the window

% load calibration
load(calibration_file_name)
window_center=Calibration.window.center_coords+offset_correction;
offset=round(window_center*1000);

% take background image
im=double(imread(im_name));

switch 2
    case 0
    case 1
        im=flipud(im(:,:,2));
        im=imresize(im,10);
    case 2
        %%
        im_new=im;
        im_new_resized=zeros([size(im_new,1)*10 size(im_new,1)*10 3]);
        for iPlane=1:3
            im_new(:,:,iPlane)=flipud(im_new(:,:,iPlane));
            im_new_resized(:,:,iPlane)=imresize(im_new(:,:,iPlane),10);
        end
        im=im_new_resized/256;
end
figure(39)
clf

H=imagesc(im);
hold on
colormap(green)
axis xy
axis equal
R=Calibration.window.radius*1e3;
axis([offset(1)-R offset(1)+R offset(2)-R offset(2)+R])

% get one example session
dataset_folder=fullfile(dataset_root,animal_ID);

files=scandir(dataset_folder,'mat');
nFiles=length(files);

parameter='EL';
direction=1;
switch parameter
    case 'AZ'
        parameter_index=3;
        nConditions=8;
        direction=-1;
    case 'EL'
        parameter_index=4;
        nConditions=4;
    case 'size'
        parameter_index=5;
        nConditions=15;
end

color_scheme=jet(nConditions);
if direction==-1
    color_scheme=flipud(color_scheme);
end

%%
for iFile=1:nFiles
    load_name=fullfile(dataset_folder,files(iFile).name);
    
    load(load_name,'dataset');
    D(iFile)=dataset;
    
    % get data from dataset
    FOV_coord=offset+dataset.FOV_info.center-[-.4 .5]*1e3;
    FOV_size=dataset.FOV_info.size_um;
    MIP=dataset.MIP_std.data;
    
    % compute rect
    FOV_rect=[0 0 FOV_size];
    %ROI=CenterRectOnPoint(FOV_rect,FOV_coord(1),FOV_coord(2));
    
    % prepare MIP
    MIP_scaled=flipud(imresize(MIP,[FOV_rect(4) FOV_rect(3)]));
    MIP_scaled=calc_gamma(MIP_scaled,.5);
    MIP_scaled=MIP_scaled-min(MIP_scaled(:));
    MIP_scaled=MIP_scaled/max(MIP_scaled(:));
    
    % place MIP in BG im
    lower_left=floor([FOV_coord(1)-FOV_rect(3)/2 FOV_coord(2)-FOV_rect(4)/2]);
    coords=round([lower_left lower_left]+FOV_rect);
    patch=im(coords(2)+1:coords(4),coords(1)+1:coords(3));
    MIP_scaled=MIP_scaled*max(patch(:));
    %im(coords(2)+1:coords(4),coords(1)+1:coords(3))=MIP_scaled;
    
    % get ROI definitions
    ROIs=dataset.ROI_definitions;
    nROIs=length(ROIs);
    ROI_coords=cat(1,ROIs.center_coords);
    
    %% For each ROI, do RF analysis and plot either azimuth or elevation
    % colored ROI areas
    ROI_vector=cat(1,ROIs.ROI_nr);
    TH=.4;
    dataMatrix=zeros(nROIs,5);
    for iROI=1:nROIs
        % perform analysis per cell
        ROI=ROIs(iROI);
        ROI_nr=ROI_vector(iROI);
        results=dataset.RF_analysis(iROI,0);
        map=results.condAverage.RF_map;
        map_TH=map>TH;
        
        % extract value of parameter of interest
        if any(map_TH(:)==1)
            Regions=regionprops(map_TH,'Area','Centroid');
            Region_area=cat(1,Regions.Area);
            if any(Region_area>1)
                sel=find(Region_area==max(Region_area),1,'first');
                RF=Regions(sel);
                AZ=RF.Centroid(1);
                EL=RF.Centroid(2);
                RF_size=RF.Area;
                if RF_size>nConditions
                    RF_size=nConditions;
                end
            else
                AZ=NaN;
                EL=NaN;
                RF_size=NaN;
            end
        else            
            AZ=NaN;
            EL=NaN;
            RF_size=NaN;
        end
        dataMatrix(iROI,:)=[iROI ROI_nr AZ EL RF_size];
        
        value=dataMatrix(iROI,parameter_index);
        
        % show on background image
        mask=poly2mask(ROI.coords_MIP(:,1),ROI.coords_MIP(:,2),dataset.FOV_info.size_px(2),dataset.FOV_info.size_px(1));
        mask=flipud(imresize(mask,[FOV_rect(4) FOV_rect(3)]));        
        patch=im(coords(2)+1:coords(4),coords(1)+1:coords(3),:);
        
        if isnan(value)
            color_vector=[0 0 0];
        else
            % pick a color
            color_vector=color_scheme(floor(value),:);
            
            mask=imdilate(mask,strel('disk',15));
            for iPlane=1:3
                patch_temp=im(coords(2)+1:coords(4),coords(1)+1:coords(3),iPlane);                                
                patch_temp(mask==1)=color_vector(iPlane);
                patch(:,:,iPlane)=patch_temp;
            end
            
            im(coords(2)+1:coords(4),coords(1)+1:coords(3),:)=patch;
        end
    end
    %dataMatrix
    
    
    
    
    %     %% get ROI shaped mask
    %     nROIs=length(ROIs);
    %     for iROI=1:nROIs
    %         ROI=ROIs(iROI);
    %         mask=poly2mask(ROI.coords_MIP(:,1),ROI.coords_MIP(:,2),dataset.FOV_info.size_px(2),dataset.FOV_info.size_px(1));
    %         mask=flipud(imresize(mask,[FOV_rect(4) FOV_rect(3)]));
    %
    %         patch=im(coords(2)+1:coords(4),coords(1)+1:coords(3));
    %         %patch(mask==1)=patch(mask==1)+25;
    %         patch(mask==1)=255;
    %         im(coords(2)+1:coords(4),coords(1)+1:coords(3))=patch;
    %     end
    
    
    %%
    %plot(FOV_coord(1),FOV_coord(2),'m*')
    %plotRect(ROI,'k');
    %plot(lower_left(1),lower_left(2),'co')
    %plot(coords(3),coords(4),'ys')
    set(H,'cData',im)
    %set(gcf,'cLim','on')
    %set(gca,'cLim',[0 2000])
    %plot(lower_left(1)+ROI_coords(:,1)*dataset.FOV_info.pixel_size_micron(2),lower_left(2)+(FOV_rect(4)-ROI_coords(:,2)*dataset.FOV_info.pixel_size_micron(1)),'r.')
    
    drawnow
end
%%% 

%%
figure(45)
switch parameter
    case 'AZ'
        color_values=repmat(1:nConditions,4,1);
        label_str='AZIMUTH';
    case 'EL'
        color_values=repmat(1:nConditions,8,1)';
        label_str='ELEVATION';
    case 'size'
end
imagesc(color_values)
colormap(color_scheme)
axis xy
xlabel(label_str)

%%



if 0
    %%
    im_folder=fullfile(dataset_folder,'renders')
    save_name=exp_name;
    save_name=fullfile(im_folder,strrep(save_name,'/','_'));
    savec(save_name)
    im_save=real(calc_gamma(im,.7));
    im_save=im_save/max(im_save(:))*256;
    
    im_save=uint8(im_save);
    if size(im_save,3)==3
        for iPlane=1:3
           im_save(:,:,iPlane)=flipud(im_save(:,:,iPlane)); 
        end
    else
        im_save=flipud(im_save);
        im_save=cat(3,im_save*0,im_save,im_save*0); % write to the green channel
    end
    
    %%% crop image
    avg=mean(im_save,3);
    col_avg=mean(avg,1);
    col_range=[find(abs(diff(col_avg>0))==1,1,'first') find(abs(diff(col_avg>0))==1,1,'last')];

    row_avg=mean(avg,2);
    row_range=[find(abs(diff(row_avg>0))==1,1,'first') find(abs(diff(row_avg>0))==1,1,'last')];
    
    im_crop=im_save(row_range(1):row_range(2),col_range(1):col_range(2),:);
    %size(im_crop)
    %figure(3);
    %plot(col_avg>0)
    %imshow(im_crop,[])
    
    %%% save image to file
    %die
    imwrite(im_crop,[save_name '_' parameter '_onlySelective_largeROI.png'])
end


%% put FOV in place
%plot(FOV_coord(1),FOV_coord(2),'m*')



% plot ROIs overlay
