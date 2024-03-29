%clear all
%clc

%%
switch 5
    case 1
        data_folder='D:\Dropbox (coxlab)\2p-data\2015-03-05_AF11\resaved';
        session_nr=18;
        file_name=fullfile(data_folder,sprintf('20150305_AF11_%03d_resaved.tif',session_nr));
        %file_name=fullfile(data_folder,sprintf('20150305_AF11_%03d.tif',session_nr));
    case 2
        session_nr=1;
        data_folder='D:\Dropbox (coxlab)\2p-data\2015-08-03_AH02_init\resaved';
        file_name=fullfile(data_folder,sprintf('2015-08-03_AH02_%03d_resaved.tif',session_nr))
    case 3
        session_nr=1;
        data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-03_AH02_init/resaved';
        file_name=fullfile(data_folder,sprintf('2015-08-03_AH02_%03d_resaved.tif',session_nr))
    case 4
        session_nr=1;
        data_folder='C:\Users\labuser\Documents\ben\Movies\2015-08-06_test';
        file_name=fullfile(data_folder,sprintf('test_%03d.tif',session_nr))
    case 5
        session_nr=28;
        data_folder='C:\Users\labuser\Documents\ben\Movies\2015-08-10_AH02';
        file_name=fullfile(data_folder,sprintf('2015-08-10_AH02_%03d.tif',session_nr))
end
%file_name='C:\Users\labuser\Desktop\20150305_AF11_009_resaved.tif';



info=imfinfo(file_name);

nFrames=length(info);
W=info(1).Width;
H=info(1).Height;
dataMatrix=zeros(nFrames,7);

%%
offset=double(intmax('uint16')/2);
for iFrame=5%1:nFrames %min([500 nFrames])
    % read out last line of the frame, here we store frame specific
    % info: bitCodes, position, laser power
    frame=imread(file_name,iFrame,'info',info);
    figure(66)
    imshow(frame,[])
    
    flyback_line=double(imread(file_name,iFrame,'info',info,'PixelRegion',{[H H],[1 W]}));
    parse_flyback_line(flyback_line)
end
%%
%plot(flyback_line)