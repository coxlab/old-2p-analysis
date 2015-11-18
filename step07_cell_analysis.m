clear all
clc

header_script

animal_ID='AH03';

dataset_folder=fullfile(dataset_root,animal_ID,'cell_data_files');
dataset_files=scandir(dataset_folder,'mat');
nFiles=length(dataset_files);

offset_correction=[0 0];
switch animal_ID
    case 'AH02'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH02_20150803.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH02_resaved_im.png';
    case 'AH03'
        if ispc
            calibration_file_name='C:\Users\LBP\Documents\GitHub\MotionGUI\Calibrations\AH03_20150807.mat';
            im_name='C:\Users\LBP\Documents\GitHub\MotionGUI\Images\2015-08-10_AH03_im.png';
        else
            calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH03_20150807.mat';
            im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH03_im.png';
        end
        %offset_correction=[-.4 .5];
    case 'AH05'
        if ispc
            calibration_file_name='C:\Users\LBP\Documents\GitHub\MotionGUI\Calibrations\AH05_20150901.mat';
            im_name='C:\Users\LBP\Documents\GitHub\MotionGUI\Images\2015-09-01_AH05_im.png';
        else
            calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH05_20150901.mat';
            im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-09-01_AH05_im.png';
        end
end



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
resize_factor=.1;
BG=double(imread(im_name))/256;
%BG=imresize(BG,resize_factor);
BG=flipud(BG);
figure(333)
clf
imshow(BG,[])
colormap(green)
hold on
marker_size=8;
p(1)=plot(0,0,'k.','markerSize',marker_size);
p(2)=plot(0,0,'.','color',[1 1 1]*.4,'markerSize',marker_size);
p(3)=plot(0,0,'r.','markerSize',marker_size);
t=title('Variable name');
axis equal
axis xy
axis([400 1000 400 1000])
drawnow



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
responsive_positions=cat(1,cell_data.nResponsive_positions);
%sel=responsive_cells>0;
%[sum(sel) length(sel)]

%size(cat(1,cell_data.RF_center))

RF_center=cat(1,cell_data.RF_center);
AZ=RF_center(:,1);
EL=RF_center(:,2);
RF_sizes=cat(1,cell_data.RF_size);
cell_size=cat(1,cell_data.cell_size);

%
sparseness_avg=cat(1,cell_data.sparseness_avg);
sparseness_max=NaN(size(sparseness_avg));
for iCell=1:nCells
    if cell_data(iCell).nResponsive_positions>0
        sparseness_max(iCell)=max(cell_data(iCell).sparseness_per_position);
    end
end
invariance_avg=cat(1,cell_data.invariance_avg);
%%

analysis_variable_names={'Responsive','Selective positions','Azimuth','Elevation','RF size','Cell Size','Sparseness','Invariance'};
analysis_variable=4;

%sel1=false(nCells,1);
sel1=responsive_positions>0;
%sel1=responsive_cells>0&cell_size>10*10;
switch analysis_variable
    case 1
        sel2=responsive_positions>0;
        %sel2=responsive_positions>0;
    case 2        
        sel2=responsive_positions>4;
        %figure(334)
        %MAP_avg=mean(cat(3,cell_data(sel2).RF_map_TH),3);
        %imagesc(flipud(MAP_avg))
        %xlabel('Periphery <=> Center')
    case 3
        sel2=sel1==1&AZ<=4; % periphery to center 1-8
        %sel2=sel1==1&AZ>=9 - 5; % center to periphery 
    case 4
        sel2=sel1==1&EL<=2;
    case 5 % Receptive field size
        sel2=sel1==1&RF_sizes>2;
    case 6 % somas and not dendrites
        sel2=sel1==1&cell_size>30*10;
    case 7 % selective cells
        sel2=sel1==1&sparseness_avg>3/10;
    case 8 % tolerant cells
        sel2=sel1==1&invariance_avg>2/10;
        
    case 9
        %V1_coords=[5899 8957 ; 5888 6444 ; 7521 8957];
        V1_coords=[5790 6160 ; 5790 9140 ; 7770 9140];
        sel2=inpolygon(cell_locations(:,1),cell_locations(:,2),V1_coords(:,1),V1_coords(:,2));
        
        % combos
    case 10 
        %sel1=responsive_cells>0&cell_size>20*10;
        sel2=sel1==1&AZ<=1;
    case 11 
        %sel1=responsive_cells>0&cell_size>12*10;
        sel2=sel1==1&AZ>=5;
    case 12
        %sel1=responsive_cells>0&cell_size>12*10;
        sel2=sel1==1&EL>=4;        
end

[length(sel1) sum(sel1) sum(sel2)]
tabulate(sel1)
tabulate(sel2(sel1))

%%%
%figure(333)
%plot(cell_data(1).offset(1)*resize_factor,cell_data(1).offset(2)*resize_factor,'rs')

%plot(FOV_centers(:,1)*resize_factor,FOV_centers(:,2)*resize_factor,'yo')
%plot(upper_left(:,1)*resize_factor,upper_left(:,2)*resize_factor,'rs')
%plot(cell_locations(:,1)*resize_factor,cell_locations(:,2)*resize_factor,'ks')

%plot(cell_locations(sel1,1)*resize_factor,cell_locations(sel1,2)*resize_factor,'k.')
%plot(cell_locations(sel2,1)*resize_factor,cell_locations(sel2,2)*resize_factor,'r.')


set(p(2),'xData',cell_locations(~sel1,1)*resize_factor,'yData',cell_locations(~sel1,2)*resize_factor)
set(p(1),'xData',cell_locations(sel1,1)*resize_factor,'yData',cell_locations(sel1,2)*resize_factor)
set(p(3),'xData',cell_locations(sel2,1)*resize_factor,'yData',cell_locations(sel2,2)*resize_factor)
set(t,'string',analysis_variable_names{analysis_variable})
axis([570 920 590 920])
%hold off

%% plot the percentage of neuron within a certain area, relative to scambled values
% to do

% decide to work in image space or physical space

which_space=2;
switch which_space
    case 1
        % create a big empty matrix
        
        % fill pixel at each cell_location with value of that cell
        
        % smooth this jittered map, ignore NaNs in between
        
        % try interp smoothed surface
        X=round(cell_locations(sel2,1));
        Y=round(cell_locations(sel2,2));
        Z=AZ(sel2);
        plot(X,Y,'.')
        
        S=sparse(X-min(X)+1,Y-min(Y)+1,Z);
        spy(S)
        
        F=full(S);
        F_s=imresize(F,.2);
                        
        imshow(rot90(F_s),[])
        colormap jet
    case 2
        %% Set parameters
        resample_factor=1;
        stride_length=10*resample_factor;
        
        selection_region='circle';
        window_size=250; % window size in which to look for cells if selection_region is rect
        inclusion_radius=125; % if selection_region is circle
        nCells_min=3; % minimal number of cells before doing analysis
                
        X=cell_locations(:,1);
        Y=cell_locations(:,2);
        range_x=[min(X) max(X)];
        range_y=[min(Y) max(Y)];
        nSteps_x=round(diff(range_x)/stride_length);
        nSteps_y=round(diff(range_y)/stride_length);
        [rows,cols]=meshgrid(linspace(range_x(1),range_x(2),nSteps_x),linspace(range_y(1),range_y(2),nSteps_y));
        
        R=round(rows(:));
        C=round(cols(:));
        
        N=length(R);        
        res_image_vector=zeros(N,1);
        pos_rect=[0 0 window_size window_size];
        
        if analysis_variable==1
            pre_selection=responsive_positions>-1;
        else
            pre_selection=responsive_positions>0;
        end
        nSelected_neurons=sum(pre_selection);
        
        tic
        for iPos=1:N
            switch selection_region
                case 'rect'
                    rect=CenterRectOnPoint(pos_rect,R(iPos),C(iPos));
                    sel=inpolygon(X,Y,rect([1 3]),rect([2 4]))&pre_selection;
                case 'circle'
                    [theta,rho]=cart2pol(X-R(iPos),Y-C(iPos));
                    sel=rho<inclusion_radius&pre_selection;
            end            
            if sum(sel)>nCells_min                
                switch analysis_variable
                    case 1 % %responsive                        
                        res_image_vector(iPos)=mean(responsive_positions(sel)>0)*100;
                    case 2 % #responsive positions
                        res_image_vector(iPos)=mean(responsive_positions(sel));                        
                    case 3 % Azimuth
                        res_image_vector(iPos)=mean(9-AZ(sel));
                    case 4 % elevation
                        res_image_vector(iPos)=mean(4-EL(sel));
                    case 5 % receptive field size
                        res_image_vector(iPos)=mean(RF_sizes(sel));
                    case 6 % sparseness                         
                        res_image_vector(iPos)=mean(sparseness_max(sel))*100;
                    case 7 % invariance
                        res_image_vector(iPos)=mean(invariance_avg(sel))*100;                        
                end
            end
        end    
        toc
        
        %% Show map
        res_image=reshape(res_image_vector,nSteps_y,nSteps_x);
        res_image_smooth=imgaussfilt(res_image,3);
        
        figure(456)
        imagesc((res_image_smooth))
        axis xy
        axis square
        set(gca,'XTickLabel',get(gca,'XTick')*stride_length)
        set(gca,'YTickLabel',get(gca,'YTick')*stride_length)
        xlabel('Medial-Lateral position (µm)')
        ylabel('Posterior-Anterior position (µm)')
        title(sprintf([analysis_variable_names{analysis_variable} ' (%d out of %d neurons)'],[nSelected_neurons nCells]))
        set(gca,'CLim',[min(res_image(res_image(:)>0)) max(res_image(res_image(:)>0))])
        colormap(jet)
        colorbar        
        
        if 0
            %%
            file_name='/Users/benvermaercke/Desktop/test.png';
            
            A=uint8(res_image_smooth/max(res_image_smooth(:))*256);
            
            
            imwrite(im,file_name)
        end
end




%% Divide into V1 and non-V1, using in poly
%V1_coords=[5899 8957 ; 5888 6444 ; 7521 8957];
V1_coords=[5790 6160 ; 5790 9140 ; 7770 9140];
V1_selection=inpolygon(cell_locations(:,1),cell_locations(:,2),V1_coords(:,1),V1_coords(:,2));

cell_selection=responsive_positions>1&cell_size>0;

%% Receptive field sizes 
% M=[responsive_positions V1_selection RF_sizes];
% M=M(responsive_positions>0,:);

% disp('receptive field size')
% p_value=kruskalwallis(M(:,3),1-M(:,2))

M=[responsive_positions 1-V1_selection RF_sizes];
M=M(cell_selection,:);
N=size(M,1);

m_V1=mean(M(M(:,2)==0,3));
m_LM=mean(M(M(:,2)==1,3));
p_value=kruskalwallis(M(:,3),M(:,2),'off');
fprintf('RF_sizes (V1:%3.2f vs. LM:%3.2f) p=%3.2f (N=%d) \n',[m_V1 m_LM p_value N])


% p=0.5077

%% Sparseness, proxy for selectivity

% V1
% sel1=responsive_positions>0&V1_selection==1;
% sel2=sel1==1&sparseness_avg>3/10;
% tabulate(sel2(sel1))
% 
% % extrastriate
% sel1=responsive_positions>0&V1_selection==0;
% sel2=sel1==1&sparseness_avg>3/10;
% tabulate(sel2(sel1))

% plot
% subplot(211)
% hist(sparseness_avg(responsive_positions>0&V1_selection==1),50)
% axis([0 1 0 100])
% subplot(212)
% hist(sparseness_avg(responsive_positions>0&V1_selection==0),50)
% axis([0 1 0 100])


% stats
M=[responsive_positions 1-V1_selection sparseness_avg];
M=M(cell_selection,:);
N=size(M,1);

m_V1=mean(M(M(:,2)==0,3));
m_LM=mean(M(M(:,2)==1,3));
% tabulate(M(M(:,2)==0,3)>.3)
% tabulate(M(M(:,2)==1,3)>.3)

p_value=kruskalwallis(M(:,3),M(:,2),'off');
fprintf('Sparseness avg (V1:%3.2f vs. LM:%3.2f) p=%3.2f (N=%d) \n',[m_V1 m_LM p_value N])


% p=0.2804, no difference in sparseness!

% stats sparseness_max
M=[responsive_positions 1-V1_selection sparseness_max];
M=M(cell_selection,:);
N=size(M,1);

m_V1=mean(M(M(:,2)==0,3));
m_LM=mean(M(M(:,2)==1,3));
% tabulate(M(M(:,2)==0,3)>.3)
% tabulate(M(M(:,2)==1,3)>.3)

p_value=kruskalwallis(M(:,3),M(:,2),'off');
fprintf('Sparseness Max (V1:%3.2f vs. LM:%3.2f) p=%3.2f (N=%d) \n',[m_V1 m_LM p_value N])




%% invariance
M=[responsive_positions 1-V1_selection invariance_avg];
M=M(cell_selection,:);
N=size(M,1);

m_V1=mean(M(M(:,2)==0,3));
m_LM=mean(M(M(:,2)==1,3));
p_value=kruskalwallis(M(:,3),M(:,2),'off');
fprintf('Invariance (V1:%3.2f vs. LM:%3.2f) p=%3.2f (N=%d) \n',[m_V1 m_LM p_value N])



%%% Story: 
% for cells with at least 2 RF positions
% median RF size is slightly bigger (ns. p=0.51)
% median sparseness at best position is lower (p=0.02, ns. if averaged over positions p=0.39)
% median invariance over positions is higher (p<1e-5)


