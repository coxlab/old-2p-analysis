clear all
clc

header_script

write_MIPs=1;
write_substacks=0;
start_index=1;
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
    case '2015-08-26_AH06'
        iFile=2;
        pixel_size_micron=[500 680]./[191 512];
    case '2015-09-01_AJ01'
        switch 2
            case 1
                iFile=3;
                pixel_size_micron=[500 680]./[191 512];
            case 2 % in case we want to stitch 1 and 2 using FIJI
                iFile=2;start_index=29;
                pixel_size_micron=[500 680]./[191 512];
        end
        
    case '2015-09-01_AH05'
        iFile=1;
        pixel_size_micron=[500 680]./[191 512];
        
    case '082715 SR101 vessels imaging'
        %iFile=1;start_index=1; session_vector=1:21;
        %iFile=2;start_index=21; session_vector=2:15;
        iFile=3;start_index=35; session_vector=2:22;
        
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
    
    if ~exist('session_vector','var')
        session_vector=1:nTrajectories;
    end
    
    %     plot(zscore(z_depth))
    %     hold on
    %     plot(zscore(diff(z_depth)))
    %     plot(zscore(A))
    %     hold off
    %%
    xy=zeros(nTrajectories,2);
    for iTrack=1:nTrajectories
        if ismember(iTrack,session_vector)
            sel=trajectory_parts==iTrack;
            idx=find(sel);
            xy(iTrack,:)=M(idx(1),1:2)/1000;
            
            
            save_folder=fullfile(session_data.folder_info.save_folder,'substacks',sprintf('session%02d',iFile));
            frames=session_data.get_frames(idx);
                
            %idx([1 end])'
            if write_substacks==1                                
                %%% write to tif stack                
                tif_name=fullfile(save_folder,sprintf('substack_%03d.tif',start_index-1+iTrack));
                savec(tif_name)
                session_data.export_movie(tif_name,frames)
            end
            
            if write_MIPs==1
                %%% save avg projection
                frames_avg=imresize(mean(frames,3),[191*2 512]);
                tif_name=fullfile(save_folder,sprintf('MIP_%03d.tif',start_index-1+iTrack));
                savec(tif_name)
                session_data.export_movie(tif_name,frames_avg)
            end
        else
            fprintf('Skipping track %d\n',iTrack)
        end
    end
    
    %%
    if ismac
        %session_data.imshow(frames)
        plot(xy(:,1),xy(:,2),'o-')
        %start_index=35
        [(1:size(xy,1))' (1:size(xy,1))'-1+start_index xy]
    else
        [(1:size(xy,1))' xy]
    end
    
    %[max_vals,min_vals]=localMaxMin(z_depth)
end

%%% Run create_backgroundImage_from_volumeStack.m after this 

