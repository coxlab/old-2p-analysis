clear all
clc

header_script
animal_ID='AH03';
plot_it=0;
save_data=1;

TH=2;

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
        calibration_file_name='/Users/benvermaercke/CoxLab/MotionGUI/Calibrations/AH05_20150901.mat';
        im_name='/Users/benvermaercke/CoxLab/MotionGUI/Images/2015-09-01_AH05_im.png';
end

dataset_folder=fullfile(dataset_root,animal_ID);
dataset_files=scandir(dataset_folder,'mat');
nFiles=length(dataset_files);


switch plot_it
    case 0
    case 1
    case 2
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
    otherwise
end

%%
t0=clock;
for iFile=1:nFiles
    load_name=fullfile(dataset_folder,dataset_files(iFile).name);
    save_name=fullfile(dataset_folder,'cell_data_files',['cell_data_' dataset_files(iFile).name]);
    if exist(load_name,'file')
        load(load_name,'dataset')
        
        var=4;
        ROI_definitions=dataset.ROI_definitions;
        nROIs=length(ROI_definitions);
        for iROI=1:nROIs
            
            %%% Create a cell object based on the cell_data class, to hold
            %%% all data and perform all subsequent analyses.
            cell_data(iROI)=cell_processor(iROI,dataset);
            
            cell_data(iROI).set_coordinate_frame(calibration_file_name,offset_correction)
            cell_data(iROI).get_cell_location()
            
            %%% Build condition matrix
            cell_data(iROI).build_condition_matrix(dataset.STIM)
            %trace=dataset.RESP(:,iROI);
            trace=dataset.SPIKE(:,iROI);
            cell_data(iROI).add_trace(trace)
            
            
            cell_data(iROI).do_RF_analysis()
            cell_data(iROI).do_RF_analysis_shuffled()
            %%
            cell_data(iROI).do_threshold(TH)
            
            %cell_data(iROI).show_RF_map(cell_data(iROI).RF_map_TH)
            
            %%
            cell_data(iROI).do_stimSelect_analysis()
            cell_data(iROI).get_sparseness()
            cell_data(iROI).calc_invariance()
            
            %% special methods
            % find full map for most responsive stimulus, or any stim if
            % argument is provided
            cell_data(iROI).calc_invariance_single_stimulus()
        end
                        
        %% evaluation
        switch plot_it
            case 0
            case 1
                %%
                subplot(121)
                INV=cat(1,cell_data.invariance_avg);
                hist(INV)
                
                %%
                subplot(122)
                RF_all=cat(3,cell_data.RF_map_TH);
                imagesc(mean(RF_all,3))
                colormap parula
                axis equal
                axis tight
                axis xy
            case 2
                window_center=cell_data(1).offset*resize_factor;
                FOV_center_abs=(window_center+cell_data(1).FOV_info.center*resize_factor) - offset_correction*100;
                plot(window_center(1),window_center(2),'wo')
                plot(FOV_center_abs(1),FOV_center_abs(2),'ms')
                aperture=2.5e3*resize_factor;
                set(gca,'Xlim',[window_center(1)-aperture window_center(2)+aperture],'Ylim',[window_center(1)-aperture window_center(2)+aperture])
                drawnow
        end
        
        %%
        if save_data==1 %% overwrite without warning
            save(save_name,'cell_data')
        end
        clear cell_data
    end
    progress(iFile,nFiles,t0)
end
