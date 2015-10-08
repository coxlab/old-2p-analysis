clear all
clc

header_script

animal_ID='AH03';

load('C:\Users\LBP\Dropbox (coxlab)\2p-datasets\AH03\cell_data_files\cell_data_2015-08-19_AH03_FOV06.mat','cell_data')
nCells=length(cell_data);

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
cell_locations=zeros(nCells,2);
for iCell=1:nCells
    cell_locations(iCell,:)=cell_data(iCell).get_cell_location();
end

%%
plot(cell_locations(:,1)*resize_factor,cell_locations(:,2)*resize_factor,'.')
%axis([-350 350 -250 250])
axis xy


hold off