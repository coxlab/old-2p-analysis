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

parameter='AZ';
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

color_scheme=jet(11);

%%
cell_data=struct('file_nr',[],'ROI_nr',[],'sparseness',[],'invariance',[]);
cell_counter=0;
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
    
    % prepare MIP
    MIP_scaled=flipud(imresize(MIP,[FOV_rect(4) FOV_rect(3)]));
    MIP_scaled=calc_gamma(MIP_scaled,.5);
    MIP_scaled=MIP_scaled-min(MIP_scaled(:));
    MIP_scaled=MIP_scaled/max(MIP_scaled(:));
    
    % place MIP in BG im
    lower_left=floor([FOV_coord(1)-FOV_rect(3)/2 FOV_coord(2)-FOV_rect(4)/2]);
    upper_left=floor([FOV_coord(1)-FOV_rect(3)/2 FOV_coord(2)+FOV_rect(4)/2]);
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
                
                %%                
                nSelectedPositions=2;
                [sorted,order]=sort(map(:),'descend');
                results_stim_sel=dataset.Analyze_stimulus_selectivity(iROI,order(1:nSelectedPositions));
                X1=cat(1,results_stim_sel.resp.iCondition);
                A=cat(1,results_stim_sel.resp.condition_mean);
                
                results_stim_sel=dataset.Analyze_stimulus_selectivity(iROI,order(nSelectedPositions+1));
                X2=cat(1,results_stim_sel.resp.iCondition);
                B=cat(1,results_stim_sel.resp.condition_mean);
                
                % sparseness
                N=length(A);
                A(A<0)=0; % rectify to avoid out of bound values
                B(B<0)=0; % rectify to avoid out of bound values
                a=( (sum(A)/N)^2 ) / ( sum((A.^2)/N ) );
                S=(1-a)/(1-1/N);
                
                % invariance
                common=intersect(X1,X2);
                A_common=A(ismember(X1,common));
                B_common=B(ismember(X2,common));
                angular_similarity=1 - acos( sum(A_common.*B_common) / (norm(A_common)*norm(B_common)) )/pi*2;
                if isnan(angular_similarity)
                    angular_similarity=0;
                end
                r=corr([A(ismember(X1,common)) B(ismember(X2,common))]);
                
                value=S*10;
                
                cell_counter=cell_counter+1;
                cell_data(cell_counter).file_nr=iFile;
                cell_data(cell_counter).ROI_nr=iROI;
                cell_data(cell_counter).abs_ROI_location=upper_left-ROI.center_coords;
                cell_data(cell_counter).selected_positions=order(1:2);
                cell_data(cell_counter).sparseness=S; % Rust and DiCarlo 2012 => Vinje and Gallant 2002
                cell_data(cell_counter).correlation=r;
                cell_data(cell_counter).angular_similarity=angular_similarity; % Rust and DiCarlo 2012
                
                
                % show on background image
                mask=poly2mask(ROI.coords_MIP(:,1),ROI.coords_MIP(:,2),dataset.FOV_info.size_px(2),dataset.FOV_info.size_px(1));
                mask=flipud(imresize(mask,[FOV_rect(4) FOV_rect(3)]));
                patch=im(coords(2)+1:coords(4),coords(1)+1:coords(3),:);
                
                color_vector=color_scheme(floor(value)+1,:);
                
                mask=imdilate(mask,strel('disk',15));
                for iPlane=1:3
                    patch_temp=im(coords(2)+1:coords(4),coords(1)+1:coords(3),iPlane);                    
                    patch_temp(mask==1)=color_vector(iPlane);
                    patch(:,:,iPlane)=patch_temp;
                end
                
                im(coords(2)+1:coords(4),coords(1)+1:coords(3),:)=patch;
                set(H,'cData',im)
                
                %c=cell_data(cell_counter).abs_ROI_location;
                %plot(c(1),c(2),'o','markerSize',10,'color',color_scheme(floor(angular_similarity*10)+1,:))
                
                plot_it=0;
                if plot_it==1
                    subplot(121)
                    bar(X1,A)
                    title(S)
                    subplot(122)
                    bar(X2,B)
                    title(r(1,2))
                    drawnow
                end
            end
        end
    end        
    
    drawnow
end


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
    parameter='sparseness';
    imwrite(im_crop,[save_name '_' parameter '_onlySelective_2pos.png'])
end

%% plot overview stats
%%% define V1 vs rest
coords=cat(1,cell_data.abs_ROI_location);
V1_coords=[5899 8957 ; 5888 6444 ; 7521 8957];
V1_selection=inpolygon(coords(:,1),coords(:,2),V1_coords(:,1),V1_coords(:,2));

TE_coords=[7464 6033 ; 8572 8591 ; 8823 6044];
TE_selection=inpolygon(coords(:,1),coords(:,2),TE_coords(:,1),TE_coords(:,2));


LM_selection=~V1_selection&~TE_selection;

sum([V1_selection LM_selection TE_selection])

%plot(coords(V1_selection,1),coords(V1_selection,2),'c+')
%plot(coords(TE_selection,1),coords(TE_selection,2),'r+')



%%
A=cat(1,cell_data.sparseness);
B=cat(1,cell_data.angular_similarity);

%%% 
figure(1)
subplot(211)
bar([mean(A(V1_selection)) mean(A(LM_selection)) mean(A(TE_selection))])
subplot(212)
bar([mean(B(V1_selection)) mean(B(LM_selection)) mean(B(TE_selection))])


%%% Plot correction sparseness tolerance
figure(2)
plot(A,B,'.')
hold on
plot(A(V1_selection),B(V1_selection),'r.')
plot(A(LM_selection),B(LM_selection),'b.')
plot(A(TE_selection),B(TE_selection),'g.')
hold off
axis equal
box off
axis square
title(sprintf('Correlation sparseness and tolerance: %3.2f',corr(A,B)))

corr(A(V1_selection),B(V1_selection))
corr(A(LM_selection),B(LM_selection))
corr(A(TE_selection),B(TE_selection))


