clear all
clc


header_script

files=scandir(data_folder,'tif');
nFiles=length(files);
for iFile=9
    file_name=fullfile(data_folder,files(iFile).name);
    if exist(file_name,'file')==2
        save_name=fullfile(data_folder,'data_analysis',files(iFile).name);
        save_name=strrep(save_name,'tif','mat');
        load(save_name)
        
        session_data.rebase(data_root)
        
        %%
        session_data.get_mov_info()
        session_data.get_mov_info()
        session_data.get_scim_data()
        session_data.read_flyback()
        session_data.get_FOV_info(.85)                
        session_data.save_data()
        
        %%
        M=[cat(1,session_data.frame_info.xyz_submicron) cat(1,session_data.frame_info.laser_power)];
        
        %plot(M(:,3:4))        
        z_depth=M(:,3);
        
        A=diff(z_depth)>0;
        trajectory_parts=parse_conditions(A);
        max_frames=min(trajectory_parts(:,4));
        
        nTrajectories=size(trajectory_parts,1);
        
        write_substacks=0;
        xy=zeros(nTrajectories,2);
        for iTrack=1:nTrajectories            
            xy(iTrack,:)=M(trajectory_parts(iTrack,2),1:2)/1000;            
            
            if write_substacks==1
                %%% write to tif stack
                frames=session_data.get_frames(trajectory_parts(iTrack,2):trajectory_parts(iTrack,2)-1+max_frames);
                tif_name=fullfile(session_data.folder_info.save_folder,'substacks',sprintf('substack_%03d.tif',iTrack));
                savec(tif_name)
                %session_data.export_movie(tif_name,frames)
            end
        end
        %session_data.imshow(frames)
        plot(xy(:,1),xy(:,2),'o-')
        
        %[max_vals,min_vals]=localMaxMin(z_depth)
    end
end

