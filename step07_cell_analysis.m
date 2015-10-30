clear all
clc

header_script

animal_ID='AH03';

dataset_folder=fullfile(dataset_root,animal_ID,'cell_data_files');
dataset_files=scandir(dataset_folder,'mat');
nFiles=length(dataset_files);

if ispc
    im_name='C:\Users\LBP\Documents\GitHub\MotionGUI\Images\2015-08-10_AH03_im.png';
else
    im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH03_im.png';
end

%%
resize_factor=.1;
BG=double(imread(im_name))/256;
%BG=imresize(BG,resize_factor);
BG=flipud(BG);
figure(333)
clf
imshow(BG,[])
colormap(green)
hold on
p(1)=plot(0,0,'k.');
p(2)=plot(0,0,'r.');
axis equal
axis xy
axis([400 1000 400 1000])
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

FOV_mapping=getMapping({cell_data.session_date});
cell_locations=cat(1,cell_data.cell_location_FOV_um);
responsive_cells=cat(1,cell_data.nResponsive_positions);
%sel=responsive_cells>0;
%[sum(sel) length(sel)]

%size(cat(1,cell_data.RF_center))

RF_center=cat(1,cell_data.RF_center);
AZ=RF_center(:,1);
EL=RF_center(:,2);
RF_sizes=cat(1,cell_data.RF_size);
cell_size=cat(1,cell_data.cell_size);

sparseness_avg=cat(1,cell_data.sparseness_avg);
invariance_avg=cat(1,cell_data.invariance_avg);

switch 10
    case 1
        sel=responsive_cells>0;
    case 2        
        sel=FOV_mapping==7;        
        figure(334)
        MAP_avg=mean(cat(3,cell_data(sel).RF_map_TH),3);
        imagesc(flipud(MAP_avg))
        xlabel('Periphery <=> Center')
    case 3
        sel=AZ<=2;
    case 4
        sel=EL<=2;
    case 5
        sel=RF_sizes>2;
    case 6
        sel=cell_size>30*10;
    case 7
        sel=sparseness_avg>.3;
    case 8
        sel=invariance_avg>.1;
        
    case 10 % combos
        sel1=responsive_cells>0&cell_size>8*10;
        sel2=sel1==1&AZ<=1;
end

tabulate(sel2(sel1))

%%%
%figure(333)
%plot(cell_data(1).offset(1)*resize_factor,cell_data(1).offset(2)*resize_factor,'rs')

%plot(FOV_centers(:,1)*resize_factor,FOV_centers(:,2)*resize_factor,'yo')
%plot(upper_left(:,1)*resize_factor,upper_left(:,2)*resize_factor,'rs')
%plot(cell_locations(:,1)*resize_factor,cell_locations(:,2)*resize_factor,'ks')

%plot(cell_locations(sel1,1)*resize_factor,cell_locations(sel1,2)*resize_factor,'k.')
%plot(cell_locations(sel2,1)*resize_factor,cell_locations(sel2,2)*resize_factor,'r.')


set(p(1),'xData',cell_locations(sel1,1)*resize_factor,'yData',cell_locations(sel1,2)*resize_factor)
set(p(2),'xData',cell_locations(sel2,1)*resize_factor,'yData',cell_locations(sel2,2)*resize_factor)
%axis([-350 350 -250 250])
%hold off

