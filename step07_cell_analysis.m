clear all
clc

header_script

animal_ID='AH03';

dataset_folder=fullfile(dataset_root,animal_ID,'cell_data_files');
dataset_files=scandir(dataset_folder,'mat');
nFiles=length(dataset_files);


im_name='C:\Users\LBP\Documents\GitHub\MotionGUI\Images\2015-08-10_AH03_im.png';
resize_factor=.1;
BG=double(imread(im_name))/256;
%BG=imresize(BG,resize_factor);
BG=flipud(BG);
figure(333)
clf
imshow(BG,[])
colormap(green)
hold on
axis equal
axis xy
drawnow




%%
tic
cell_data=[];
for iFile=1:nFiles    
    load_name=fullfile(dataset_folder,dataset_files(iFile).name);
    S=load(load_name);
    cell_data=cat(2,cell_data,S.cell_data);
end
toc
nCells=length(cell_data);

%%
% cell_locations=zeros(nCells,2);
% upper_left=zeros(nCells,2);
% FOV_centers=zeros(nCells,2);
% 
% for iCell=1:nCells
%     [cell_locations(iCell,:), upper_left(iCell,:), FOV_centers(iCell,:)]=cell_data(iCell).get_cell_location();    
%     
%     %%% re-threshold
%     cell_data(iCell).do_threshold(2)
% end
%%
cell_locations=cell_data.cell_location_FOV_um;
responsive_cells=cat(1,cell_data.nResponsive_positions);
RF_center=cat(1,cell_data.RF_center);
AZ=RF_center(:,1);
EL=RF_center(:,2);
RF_sizes=cat(1,cell_data.RF_size);
sel=responsive_cells>0;

tabulate(sel)

%%
%plot(cell_data(1).offset(1)*resize_factor,cell_data(1).offset(2)*resize_factor,'rs')

%plot(FOV_centers(:,1)*resize_factor,FOV_centers(:,2)*resize_factor,'yo')
%plot(upper_left(:,1)*resize_factor,upper_left(:,2)*resize_factor,'rs')
plot(cell_locations(:,1)*resize_factor,cell_locations(:,2)*resize_factor,'k.')
plot(cell_locations(sel,1)*resize_factor,cell_locations(sel,2)*resize_factor,'.')

%axis([-350 350 -250 250])
axis xy
%hold off