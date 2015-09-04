clear all
clc

%%% Optimize parameters of the Circular Hough Transform (imfindcircles)
%%% using ground truth provided by manual scoring. If not ground truth, at
%%% least we should get results that are close to our way of drawing ROIs.

%%% Define example database:
load_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-10_AH02/resaved/data_analysis/2015-08-10_AH02_001.mat';
load(load_name)


%%% Get ground truth
GT.centers=cat(1,session_data.ROI_definitions(session_data.ROI_definition_nr).ROI.center_coords);



%%% Get MIP to work with
%MIP=session_data.MIP_cc_local.data;
MIP=session_data.MIP_avg.data;
MIP_processed=MIP;

%%% Do adaptive TH
if 1
    ws=20;
    C=-.03;
    tm=0;
    MIP_processed=adaptivethreshold(MIP_processed,ws,C,tm);
end

%%% Remove speckle
SE=[0 1 0 ; 1 1 1 ; 0 1 0];
MIP_processed=imopen(MIP_processed,SE);



Sensitivity=.90; % .85
EdgeThreshold=.03; % graythresh(calc_gamma(MIP,.5))
[centers, radii]=imfindcircles(MIP_processed,[4 10],'Sensitivity',Sensitivity,'EdgeThreshold',EdgeThreshold);

%% Create mask per ROI


%% compare with ground
% number of ROIs

% location of centers


% minimal distance

if 0
    % ICP
    q=[GT.centers zeros(size(GT.centers,1),1)]'; % ref
    p=[centers zeros(size(centers,1),1)]';
    [TR, TT, ER, t]=icp(q,p,15);
    points_matched=(TR*p+repmat(TT,1,size(p,2)))';
end


% overlap between ROIs
coords=centers;
bw=bwselect(MIP_processed,coords(:,1),coords(:,2));
bw_separated=bwlabel(bw);
ROI_vector=unique(bw_separated(:));
ROI_vector=ROI_vector(ROI_vector>0);
nROI=length(ROI_vector);
window_size=[40 40];

for iROI=1:nROI
    ROI_nr=ROI_vector(iROI);
    sel=bw_separated==ROI_nr;
    neg=imerode(sel,SE);
    edge=sel-neg;
    
    [X,Y]=find(edge);
    
    ellipse_properties=fit_ellipse(Y,X);
    
    blank=session_data.blank_ROI();
    ROI=blank.ROI;
    
    ROI.ROI_nr=ROI_nr;
    ROI.base_coord=[ellipse_properties.X0_in ellipse_properties.Y0_in];
    ROI.nCoords=size(coords,1);
    ROI.coords=[Y-ROI.base_coord(1)+window_size(1)/2+1 X-ROI.base_coord(2)+window_size(2)/2+1];
    
    ROI.ellipse_properties=ellipse_properties;
    ROI.ellipse_coords=[ellipse_properties.rotated_ellipse(1,:)'-ROI.base_coord(1)+window_size(1)/2+1 ellipse_properties.rotated_ellipse(2,:)'-ROI.base_coord(2)+window_size(2)/2+1];
    ROI.ellipse_coords_centered=ROI.ellipse_coords;
    ROI.coords_MIP=ellipse_properties.rotated_ellipse';
    ROI.coords_MIP_plot=ellipse_properties.rotated_ellipse';
    
    ROI.center_coords=[ellipse_properties.X0_in ellipse_properties.Y0_in];
    
    ROI.ROI_rect=round([ROI.center_coords ROI.center_coords]+[-window_size/2+1 window_size/2]);
    
    %%% Create masks
    % Get region from image around ROI
    try
        %T=MIP(ROI.ROI_rect(2):ROI.ROI_rect(4),ROI.ROI_rect(1):ROI.ROI_rect(3));
        MIP_max=centerOnRect(MIP,size(MIP)+window_size);
        T=MIP_max(ROI.ROI_rect(2)+window_size(2)/2:ROI.ROI_rect(4)+window_size(2)/2,ROI.ROI_rect(1)+window_size(1)/2:ROI.ROI_rect(3)+window_size(2)/2);
    catch
        figure()
        session_data.imshow(sel)
        title('Failed to fit ROI')
        break
    end
    
    % Generate mask to isolate soma pixels
    %offset=ROI.ellipse_coords_centered(1,:)-window_size/2
    mask_soma=poly2mask(ROI.ellipse_coords_centered(:,1),ROI.ellipse_coords_centered(:,2),size(T,1),size(T,2));
    %figure(50)
    %ROI.ellipse_coords_centered
    %imshow(mask_soma,[])
    
    %die
    
    % Generate mask to isolate neuropil pixels
    mask_neuropil=imresize(mask_soma,sqrt(2),'nearest');
    ext=round((size(mask_neuropil)-window_size)/2);
    mask_neuropil=mask_neuropil(ext:ext+window_size(1)-1,ext:ext+window_size(2)-1);
    mask_neuropil=mask_neuropil-mask_soma;
    
    ROI.mask_soma=mask_soma;
    ROI.mask_neuropil=mask_neuropil;
    
    ROI.timeseries_soma=[];
    ROI.time_series_neuropil=[];
    
    if isfield(session_data.ROI_definitions(1).ROI,'timeseries_neuropil')
        session_data.ROI_definitions(1).ROI=rmfield(session_data.ROI_definitions(1).ROI,'timeseries_neuropil');
    end
    %%% Store ROI properties
    try
        session_data.ROI_definitions(1).ROI(iROI)=ROI;
        session_data.ROI_definitions(1).ROI(iROI)
    catch
        session_data.ROI_definitions(1).ROI(iROI)
        ROI
        error('Objects are not identical...')
    end
end
fprintf('Found %d ROIs!\n',nROI)

%% mask per ROI using roifill
for iROI=1:nROI
    coords=session_data.ROI_definitions(1).ROI(iROI).coords_MIP;
    session_data.ROI_definitions(1).ROI(iROI).mask=poly2mask(coords(:,1),coords(:,2),size(MIP,1),size(MIP,2));
end

for iROI=1:length(session_data.ROI_definitions(2).ROI)
    coords=session_data.ROI_definitions(2).ROI(iROI).coords_MIP;
    session_data.ROI_definitions(2).ROI(iROI).mask=poly2mask(coords(:,1),coords(:,2),size(MIP,1),size(MIP,2));
end
%%
session_data.plot_ROIs(1)

%% show
%session_data.imshow(MIP,.3)
%session_data.imshow(MIP_processed)
session_data.imshow(bw)
hold on
plot(GT.centers(:,1),GT.centers(:,2),'ro','markerSize',10)
plot(centers(:,1),centers(:,2),'m*')
%plot(points_matched(:,1),points_matched(:,2),'cs')

hold off



