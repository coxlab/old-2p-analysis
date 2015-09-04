clear all
clc

projection_type='MIP_cc_local';
%projection_type='MIP_std';

load('/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-07-20_AG02/data_analysis/2015-07-20_AG02_003.mat')
green=session_data.green;

IM=session_data.(projection_type).data;

figure(57)
clf
subplot(121)
imshow(IM,[])
colormap(green)


%% Threshold image
switch projection_type
    case 'MIP_std'
        ws=20;
        C=-.01;
        tm=0;
    case 'MIP_cc_local'
        ws=20;
        C=-.03;
        tm=0;
end
bw=adaptivethreshold(IM,ws,C,tm);

subplot(122)
imshow(bw,[])
colormap(green)

%% remove speckle
SE=[0 1 0 ; 1 1 1 ; 0 1 0];
bw=imopen(bw,SE);

subplot(122)
imshow(bw,[])
colormap(green)

%% get regionprops
tic
props=regionprops(bw,'Area','Centroid');
toc

%% filter spot based on area
TH.area=[100 500];
sel=between(cat(1,props.Area),TH.area);

coords=cat(1,props(sel).Centroid);
bw=bwselect(bw,coords(:,1),coords(:,2));

subplot(122)
imshow(bw,[])
colormap(green)
%%
bw_separated=bwlabel(bw);
%bw_labeled=label2rgb(bw_separated);

%% turn into ROIs
ROI_vector=unique(bw_separated(:));
ROI_vector=ROI_vector(ROI_vector>0);
nROI=length(ROI_vector);

subplot(121)
imshow(IM,[])
colormap(green)
hold on
for iROI=1:nROI
    ROI_nr=ROI_vector(iROI);
    sel=bw_separated==ROI_nr;
    neg=imerode(sel,SE);
    edge=sel-neg;
    
    [X,Y]=find(edge);
    ellipse_properties=fit_ellipse(Y,X);
        
    plot(ellipse_properties.rotated_ellipse(1,:),ellipse_properties.rotated_ellipse(2,:),'r')
    plot(coords(iROI,1),coords(iROI,2),'m*')            
end
hold off
