%%% This script extracts all the FOVs from data_analysis files and plots
%%% then on the background GCaMP image.


clear all
clc

header_script

animal_ID='AH05';

offset_correction=[0 0];
switch animal_ID
    case 'AH02'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH02_20150803.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH02_resaved_im.png';
    case 'AH03' 
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH03_20150807.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-10_AH03_im.png';
        offset_correction=[.4 -.5];
    case 'AH05'
        %calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH05_20150814.mat';
        %im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-14_AH05_im.png';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-09-01_AH05_im.png';
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH05_20150901.mat';
        
        offset_correction_conditional=[-.4 .5];
        % first session was done in window coords, with correction from
        % AH03 still active, in the second session, we reversed this
        % correction.
        % so window coords need to be shifted and coords of session 01 need
        % not to be counter corrected
        
    case 'AH06'
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH06_20150826_v2.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-08-14_AH05_im.png';
end


%%% Get bg image and calibration coords

%%% Get all sessions/datasets for this animal
% parse folder_names
sub_folders=getSubFolders(data_root);
nFolders=length(sub_folders);
count=1;
for iFolder=1:nFolders
    f=sub_folders{iFolder};
    if ~isempty(strfind(f,animal_ID))
        if isdir(fullfile(data_root,f,'data_analysis'))
            folder_list{count}=fullfile(data_root,f,'data_analysis');
            count=count+1;
        elseif isdir(fullfile(data_root,f,'resaved','data_analysis'))
            folder_list{count}=fullfile(data_root,f,'resaved','data_analysis');
            count=count+1;
        else
            fprintf('No data for this session: %s\n',f)
        end
    end
end

%% Get all session in 1 big struct S
count=1;
nFolders=length(folder_list);
for iFolder=1:nFolders
    f=folder_list{iFolder};
    session_names=scandir(f,'mat');
    nSessions=length(session_names);
    for iSession=1:nSessions
        session_name=session_names(iSession).name;
        load_name=fullfile(f,session_name);
        load(load_name,'session_data')
        if session_data.is_static_FOV()==1
            S(count)=session_data;
            count=count+1;
        else
            fprintf('skipping %s\n',session_names(iSession).name)
        end
    end
end


%% collect coords
nSessions=length(S);
coords=[];
count=1;
for iSession=1:nSessions
    [f,folder_name]=fileparts(S(iSession).folder_info.data_folder);
    folder_name
    switch folder_name        
        case '2015-08-20_AH05'
            coords(count,:)=S(iSession).FOV_info.center+offset_correction_conditional*1e3;
            count=count+1;
        otherwise
            coords(count,:)=S(iSession).FOV_info.center;
            count=count+1;
    end
    
end
%%

coords_unique=round(unique(coords,'rows','stable'))/10;

switch 2
    case 1
        %%
        S.plot_FOV
    case 2
        
        load(calibration_file_name)
        window_center=Calibration.window.center_coords;
        offset=round(window_center*100);
        offset_corrected=round((window_center-offset_correction)*100);
        
        %%        
        im=double(imread(im_name));
        im=flipud(im(:,:,2));
        session_data.imshow(im,.7)
        hold on
        plot(coords_unique(:,1)+offset(1),coords_unique(:,2)+offset(2),'rs')
        hold off
        axis xy
        axis equal
        axis([offset_corrected(1)-250 offset_corrected(1)+250 offset_corrected(2)-250 offset_corrected(2)+250])
        %axis equal
end

% Figure out location of each cell within the window

% Plot either average map over FOV location

% either individual ROI colorcoded for RF center location