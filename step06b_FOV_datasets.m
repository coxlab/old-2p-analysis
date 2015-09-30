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

files=scandir(dataset_folder,'mat');
nFiles=length(files);

for iFile=1:nFiles
    load_name=fullfile(dataset_folder,files(iFile).name);
    
    load(load_name,'dataset');
    D(iFile)=dataset;
end

%%
%{D.animal_ID}
tic
nDatasets=length(D);
coords=zeros(nDatasets,2);
FOV_sizes=coords;
MIPs=struct;
for iDataset=1:nDatasets
    coords(iDataset,:)=D(iDataset).FOV_info.center;
    FOV_sizes(iDataset,:)=D(iDataset).FOV_info.size_um;
    MIPs(iDataset).FOV=D(iDataset).MIP_std.data;
    
    deconvolve=1;
    results=D(iDataset).RF_analysis([],deconvolve);
    MIPs(iDataset).AVG_RF_map=mean(cat(3,results.condAverage.RF_map),3);
end
toc


%% apply counter correction for window coord
coords_corr=coords-repmat(offset_correction*1000,nDatasets,1);


%%
coords_explode=zeros(size(coords_corr));
nCoords=size(coords,1);
for iCoord=1:nCoords
    coords_explode(iCoord,:)=coords_corr(iCoord,:)*3;
end

%% reduce overlap between explodes
FOV_distances=squareform(pdist(coords_explode));
[F1,F2]=find(FOV_distances<500)
coords_explode

min_dist=500;
Z=linkage(coords_explode,'average');
C=cluster(Z,'cutoff',min_dist,'Criterion','distance');

%%
spread=250;
clusters=unique(C);
nClusters=length(clusters);
for iClust=1:nClusters
    sel=C==iClust;
    N=sum(sel);
    if N>1
        M=coords_explode(sel,:);        
        center_coord=mean(M,1);
        new_coords=repmat(center_coord,N,1)+[linspace(-spread,+spread*(N-1),N) ; linspace(-spread,+spread*(N-1),N)]';        
        coords_explode(sel,:)=new_coords;
    end
end
%%

clf
plot(coords(:,1),coords(:,2),'k.')
hold on
plot(coords_corr(:,1),coords_corr(:,2),'b.')
plot(coords_explode(:,1),coords_explode(:,2),'r.')
for iCoord=1:nCoords
    text(coords_explode(iCoord,1),coords_explode(iCoord,2),sprintf('#%d',C(iCoord)))
end
hold off
axis([-1 1 -1 1]*6000)



%%


load(calibration_file_name)
window_center=Calibration.window.center_coords+offset_correction;
offset=round(window_center*1000);

coords_corrected=coords_corr+repmat(offset,nDatasets,1);
coords_explode=coords_explode+repmat(offset,nDatasets,1);

%%
im=double(imread(im_name));
im=flipud(im(:,:,2));
im=imresize(im,10);
figure(39)
clf
H=imagesc(im);
axis xy
colormap(green)
axis([offset(1) offset(1) offset(2) offset(2)]+[-1 1 -1 1]*2000*2)
axis equal
drawnow

% draw points
if 0
    %%
    hold on
    plot(coords_corrected(:,1),coords_corrected(:,2),'r*')
    hold off
end

%% draw boxes
hold on
for iCoord=1:nCoords
    % Actual boxes
    center=coords_corrected(iCoord,:);
    FOV_rect=[0 0 FOV_sizes(iCoord,:)];
    ROI=CenterRectOnPoint(FOV_rect,center(1),center(2));
    plotRect(ROI,'k');
    
    % Exploded boxes
    center=coords_explode(iCoord,:);
    FOV_rect=[0 0 FOV_sizes(iCoord,:)];
    ROI=CenterRectOnPoint(FOV_rect,center(1),center(2));
    plotRect(ROI,'r');
    
    % Connecting lines
    plot([coords_corrected(iCoord,1) coords_explode(iCoord,1)],[coords_corrected(iCoord,2) coords_explode(iCoord,2)],'r-')
end
hold off
drawnow

%% Inject MIP, RF map or cells into background image
im_added=im;
for iCoord=1:nCoords
    %plot(coords(:,1),coords(:,2),'r*')
    center=coords_corrected(iCoord,:);
    FOV_rect=[0 0 FOV_sizes(iCoord,:)];
    ROI=round(CenterRectOnPoint(FOV_rect,center(2),center(1)));
    
    MIP=flipud(imresize(MIPs(iCoord).FOV,round(FOV_sizes(iCoord,[2 1]))));
    
    upper_left=floor([center(2)-size(MIP,1)/2 center(1)-size(MIP,2)/2]);
    if 0
        %%
        figure(39)
        hold on
        plot(upper_left(2),upper_left(1),'co')
        hold off
    end
    patch=im_added(upper_left(1)+1:upper_left(1)+size(MIP,1),upper_left(2)+1:upper_left(2)+size(MIP,2));
        
    im_added(upper_left(1)+1:upper_left(1)+size(MIP,1),upper_left(2)+1:upper_left(2)+size(MIP,2))=calc_gamma(MIP,.5)*40;
    
    if 0
        %%
        figure(40)
        subplot(121)
        imshow(MIP,[])
        subplot(122)
        imshow(patch,[])
        colormap(green)
    end
end

set(H,'Cdata',im_added)
drawnow

%%
%im_added_RF=im;
for iCoord=1:nCoords
    center=coords_explode(iCoord,:);
    FOV_rect=[0 0 FOV_sizes(iCoord,:)];
    ROI=round(CenterRectOnPoint(FOV_rect,center(2),center(1)));
        
    MIP=imresize(MIPs(iCoord).AVG_RF_map,round(FOV_sizes(iCoord,[2 1])));
    MIP=MIP-min(MIP(:));
    MIP=MIP/max(MIP(:))*256*.75;
    
    upper_left=floor([center(2)-size(MIP,1)/2 center(1)-size(MIP,2)/2]);   
    patch=im_added(upper_left(1)+1:upper_left(1)+size(MIP,1),upper_left(2)+1:upper_left(2)+size(MIP,2));        
    im_added(upper_left(1)+1:upper_left(1)+size(MIP,1),upper_left(2)+1:upper_left(2)+size(MIP,2))=MIP;
end

set(H,'Cdata',im_added)







