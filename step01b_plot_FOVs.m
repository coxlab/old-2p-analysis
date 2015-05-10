clear all
clc

% BV20150408: add motion correction as separate step: done
% BV20150410: might need to add another script that grabs bg_im, FOV coords
% of all session in one animal and injection coords. so we can leave it out
% of this script
header_script

%%% Manually select a folder
%cd(data_root)
%data_folder=uigetdir(data_root);
plot_it=1;
save_it=0;

%%
loadName=fullfile(data_folder,'data_analysis','session_overview.mat');
if exist(loadName,'file')==2
    load(loadName)
else
    disp('no such file')
end

%%
nSessions=length(data_sessions);

offset=double(intmax('uint16'))/2;

if strcmpi(data_folder,'/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-05_AF11')
    FOV_matrix=[1,1379,-1083,757,340,435.200000000000;2,1379,-1083,757,170,217.600000000000;3,1379,-1083,757,170,217.600000000000;4,1434,-1002,756,170,217.600000000000;5,1140,-67,768,170,217.600000000000;6,1387,-1020,776,170,217.600000000000;7,1588,-986,740,170,217.600000000000;8,1400,-1071,784,170,217.600000000000;9,1400,-1071,784,170,217.600000000000;10,1268,-1211,788,170,217.600000000000;11,1268,-1211,788,170,217.600000000000;12,1268,-1211,788,170,217.600000000000;13,1268,-1211,1041,170,217.600000000000;14,1287,-1157,777,435.200000000000,435.200000000000;15,1379,-1083,757,340,435.200000000000;16,1379,-1083,757,340,435.200000000000;17,1379,-1083,757,340,435.200000000000;18,1407,-1058,757,170,217.600000000000];
else
    FOV_matrix=zeros(nSessions,6);
    for iSess=1:nSessions
        FOV_info=data_sessions(iSess).FOV_info;
        if isfield(FOV_info,'center')
            if FOV_info.center(1)-offset<-1e4
                center=FOV_info.center;
            else
                center=FOV_info.center-offset;
            end
            FOV_matrix(iSess,:)=[iSess center FOV_info.z_depth FOV_info.size_um];
        end
    end
end


%% Get all coords in a matrix and find overlapping FOVs
FOV_matching.dist_matrix=squareform(pdist(FOV_matrix(:,2:4)));
FOV_matching.min_dist=20;

switch 2
    case 1 % assume sessions at same field were recorded back to back
        same=FOV_matching.dist_matrix<FOV_matching.min_dist;
        clusters=bwlabel(same==1,4);
    case 2 % detect same FOV across all session, preferred method
        %%
        Z=linkage(FOV_matching.dist_matrix,'average');
        C=cluster(Z,'cutoff',FOV_matching.min_dist,'Criterion','distance');
        % renumber clusters
        cluster_vector=unique(C,'stable');
        clusters=zeros(size(C));
        for iC=1:length(cluster_vector)
            sel=C==cluster_vector(iC);
            clusters(sel)=iC;
        end
end
FOV_matching.clusters=clusters;


%%
switch plot_it
    case 0
        % not plotting 
    case 1
        %%
        figure(1)
        clf
        subplot(121)
        hold on
        axis([-3 3 -3 3]*1000)
        axis equal
        circle([0 0],2.5e3,100,'k-',1);
        for iSess=1:nSessions
            row=FOV_matrix(iSess,:);
            plot(row(2),row(3),'r*')
            plotRect([row(2)-row(5)/2 row(3)-row(6)/2 row(2)+row(6)/2 row(3)+row(5)/2],'r');
            text(row(2)+row(6)/2+rand*20,row(3)+row(5)/2+rand*10,num2str(iSess));
        end
        
        subplot(122)
        %dendrogram(Z)
        imagesc(FOV_matching.dist_matrix)
        colorbar
    case 2 % will be moved to a separate script
        %% Align background image
        loadName=fullfile(data_folder,'bg_im.jpg');
        %BG=rot90(imread(loadName),2);
        BG=fliplr(imread(loadName));       
        
        warning off
        imshow(BG,[])
        axis xy
        warning off
        hold on
        p=plot([],[],'m*');
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
        
        %% Compare coords with calibration from the file
        real_diameter=4; % mm
        center=mean(coords);
        coords_center=coords-repmat(center,size(coords,1),1);
        im_radius=mean(calc_dist([coords repmat(center,size(coords,1),1)])); 
        
        %%% rescale FOV coords microns to pixels
        scaling_factor=(im_radius/1000)/(real_diameter/2);
        FOV_coords=FOV_matrix(:,2:3)*scaling_factor; % rescale to size of image
        FOV_coords=FOV_coords+repmat(center,size(FOV_coords,1),1); % move to center
        %FOV_coords(1,:)
        FOV_rect=FOV_matrix(:,[6 5])*scaling_factor; % rescale to size of image
        
        %%% rescale injection coords mm to pixels       
        show_injections=~isempty(data_folder);
        if show_injections==1
            injection_coords=fliplr([0 0 ; .43 .65 ; .48 1.21 ; 1.14 1.72 ; .45 1.97]);
            injection_root=[0 0];
            scaling_factor=(im_radius)/(real_diameter/2);
            injection_coords=injection_coords*scaling_factor;
            injection_coords=injection_coords+repmat(injection_root,size(injection_coords,1),1);
        end 
        %%% correct for incorrect calibration
        %correction_distance=[-100 150];
        correction_distance=[-.3 .6]*scaling_factor; % based on 1p vs 2p comparison
        
        %correction_distance=[-200 250];
        FOV_coords=FOV_coords+repmat(correction_distance,size(FOV_coords,1),1);
        if show_injections==1
            correction_distance=[1100 1030];
            injection_coords=injection_coords+repmat(correction_distance,size(injection_coords,1),1);
        end
        
        %%%
        warning off
        imshow(BG,[])
        warning off
        axis xy
        hold on
        plot(center(1),center(2),'r+')
        circle(center,im_radius,100,'r-',1);
        for iSess=1:nSessions
            center=FOV_coords(iSess,:);
            rect=FOV_rect(iSess,:);
            plotRect([center(1)-rect(1)/2 center(2)-rect(2)/2 center(1)+rect(1)/2 center(2)+rect(2)/2],'r');
        end
        
        if show_injections==1
            % indicate injection sites
            plot(injection_coords(:,1),injection_coords(:,2),'bv','markerfacecolor','b','markerSize',10)
        end
        hold off
        
end


%%
if save_it==1
    %%    
    save(loadName,'FOV_matching','-append')
    disp('Saved FOV data to overview file')
end

%% this could then lead to pooling of the data and/or application of the same ROI definition over different sessions
