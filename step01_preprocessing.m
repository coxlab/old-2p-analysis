clear all
clc

header_script

files=scandir(data_folder,'tif');
nFiles=length(files);
%%

for iFile=1:nFiles
    file_name=fullfile(data_folder,files(iFile).name);
    if exist(file_name,'file')==2
        fprintf('Pre-processing file %s...\n',file_name)
        save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
        save_name=strrep(save_name,'tif','mat');
        
        if exist(save_name,'file')==0
            %%% Create file based on class file
            session_data=imaging_dataset(file_name);
            session_data.save_data()
        else
            load(save_name,'session_data')
        end
        
        %%% Make sure filenames are relative to data folder on this machine
        session_data.rebase(data_root)
        
        %%% Extract info from movie
        session_data.get_mov_info()
        session_data.get_scim_data()
        session_data.read_flyback()
        %session_data.get_FOV_info(.85)
        session_data.get_FOV_info([500 680]./[191 512])
        
        
        %%% Detect invalid frames
        session_data.find_blank_frames() % due to laser power not being on
        
        session_data.save_data()
        if session_data.is_static_FOV()==0 % if recording at 1 position
            fprintf('Session does not appear to be a static FOV...\n')
        else
            %%% Extract bitcode information from two sources
            session_data.get_scim_bitCodes()
            session_data.get_MWorks_bitCodes()
            session_data.find_offset()
            if ismac
                fprintf('Offset was determined to be %d events...\n',session_data.bitCodes.offset)
            end
            
            %%% Save results so far
            session_data.save_data()
            
            %%% Motion Correction
            if use_GPU==1&&gpuDeviceCount==1
                session_data.set_smoothing_kernel_GPU()
                %session_data.reset_reference_image();
                session_data.find_reference_image_GPU()
            else
                session_data.set_smoothing_kernel()
                %session_data.reset_reference_image();
                session_data.find_reference_image()                
            end            
            
            if 0
                %%
                imshow(calc_gamma(session_data.motion_correction.reference_image.im,.5),[])
                colormap(green)
            end
            
            %session_data.reset_motion_correction();
            session_data.do_motion_correction()
            session_data.find_motion_frames(2)
            session_data.save_data()
            if 0
                %%
                plot(session_data.motion_correction.shift_matrix(:,2))
                hold on
                bar(double(session_data.motion_correction.ignore_frames))
            end
            
            if 0
                %%
                session_data.visualize_motion_correction(1)
            end
            
            %%% Calculate z-projections
            %session_data.reset_MIPs();
            session_data.do_calc_MIPs()
            if 0
                %%
                session_data.imshow(session_data.MIP_std.data)
            end
            session_data.save_data()
            
            %session_data.do_motion_detection()
            % plot(session_data.motion_proxy)
        end
    end
end