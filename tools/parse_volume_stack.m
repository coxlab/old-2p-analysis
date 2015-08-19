clear all
clc

header_script

write_substacks=1;

switch exp_name
    case '2015-08-10_AH02/resaved'
        iFile=22;
        pixel_size_micron=[529 680]./[199 512];
    case '2015-08-10_AH03'
        iFile=1;
        pixel_size_micron=[500 680]./[191 512];
    case '2015-08-14_AH05'
        iFile=2;
        pixel_size_micron=[500 680]./[191 512];
    case '2015-08-18_AH06'
        iFile=1;
        pixel_size_micron=[500 680]./[191 512];
    otherwise
        iFile=1;
        die
end


files=scandir(data_folder,'tif');
nFiles=length(files);

file_name=fullfile(data_folder,files(iFile).name);
if exist(file_name,'file')==2
    save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
    save_name=strrep(save_name,'tif','mat');
    if exist(save_name,'file')==0
        session_data=imaging_dataset(file_name);
    else
        load(save_name)
    end
    
    session_data.rebase(data_root)
    
    %%
    session_data.get_mov_info()
    session_data.get_mov_info()
    session_data.get_scim_data()
    session_data.read_flyback()
    session_data.get_FOV_info(pixel_size_micron)
    session_data.save_data()
    
    %%
    M=[cat(1,session_data.frame_info.xyz_submicron) cat(1,session_data.frame_info.laser_power)];
    
    %plot(M(:,3:4))
    z_depth=M(:,3);
    A=diff(z_depth)>0;
    %A=imopen(A,[0 1 0]);
    A=medfilt1(double(A),3);
    trajectory_parts=bwlabel(A);
    %trajectory_parts=parse_conditions(A)
    %max_frames=min(trajectory_parts(:,4));
    %nTrajectories=size(trajectory_parts,1);
    
    trajectories=unique(trajectory_parts);
    nTrajectories=max(trajectories);
    
    %     plot(zscore(z_depth))
    %     hold on
    %     plot(zscore(diff(z_depth)))
    %     plot(zscore(A))
    %     hold off
    %%
    xy=zeros(nTrajectories,2);
    for iTrack=1:nTrajectories
        sel=trajectory_parts==iTrack;
        idx=find(sel);
        xy(iTrack,:)=M(idx(1),1:2)/1000;
        
        idx([1 end])'
        if write_substacks==1
            %%% write to tif stack
            frames=session_data.get_frames(idx);
            tif_name=fullfile(session_data.folder_info.save_folder,'substacks',sprintf('substack_%03d.tif',iTrack));
            savec(tif_name)
            %session_data.export_movie(tif_name,frames)
            
            %%% save avg projection
            frames_avg=imresize(mean(frames,3),[191*2 512]);
            savec(tif_name)
            tif_name=fullfile(session_data.folder_info.save_folder,'substacks',sprintf('MIP_%03d.tif',iTrack));
            session_data.export_movie(tif_name,frames_avg)
        end
    end
    
    if ismac
        %session_data.imshow(frames)
        plot(xy(:,1),xy(:,2),'o-')
    end
    
    %[max_vals,min_vals]=localMaxMin(z_depth)
end

