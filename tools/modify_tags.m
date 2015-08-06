clear all
clc

data_folder='D:\Dropbox (coxlab)\2p-data\2015-08-03_AH02_init\resaved';
src_file_name='2015-08-03_AH02_%03d.tif';
tgt_file_name='2015-08-03_AH02_%03d_resaved.tif';

session_nr=1;
load_name=fullfile(data_folder,sprintf(src_file_name,session_nr));

src=Tiff(load_name,'r');
info=imfinfo(load_name);

info.ImageDescription

