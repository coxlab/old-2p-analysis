clear all
clc

if ismac
    root_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/';
else
    root_folder='D:\Dropbox (coxlab)\2p-data\';
end

switch 3
    case 1 % 2015-08-03_AH02_init
        data_folder=fullfile(root_folder,'2015-08-03_AH02_init/resaved');
        src_file_name='2015-08-03_AH02_%03d.tif';
        tgt_file_name='2015-08-03_AH02_%03d_resaved.tif';
        nSessions=13;
    case 2 % 2015-03-05_AF11
        data_folder=fullfile(root_folder,'2015-03-05_AF11/resaved');
        src_file_name='20150305_AF11_%03d.tif';
        tgt_file_name='20150305_AF11_%03d_resaved.tif';
        nSessions=18;
    case 3 % 2015-08-02
        data_folder=fullfile(root_folder,'2015-08-10_AH02');
        src_file_name='2015-08-10_AH02_%03d.tif';
        tgt_file_name='resaved/2015-08-10_AH02_%03d.tif';
        nSessions=26;
end


for session_nr=1:nSessions
    %%% Read correct tag
    load_name=fullfile(data_folder,sprintf(src_file_name,session_nr));
    if exist(load_name,'file')==2
        src=Tiff(load_name,'r');
        tag_content=src.getTag('ImageDescription');
        src.close()
                        
        new_tag=[tag_content(1:end-1) char(10)];
        
        if 1
            %% Write tag to target file
            load_name=fullfile(data_folder,sprintf(tgt_file_name,session_nr));
            if exist(load_name,'file')==2
                tgt=Tiff(load_name,'r+');
                tgt.setTag('ImageDescription',new_tag)
                tgt.close()
                
                %% double check
                tgt=Tiff(load_name,'r');
                tag_content=tgt.getTag('ImageDescription');
                tgt.close()
            else
                disp('No file to correct..')
            end
        end
    else
        disp('Skipping')
        load_name
    end
end
