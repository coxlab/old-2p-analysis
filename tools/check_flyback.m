clear all
clc

%%
switch 1
    case 1
        data_folder='D:\Dropbox (coxlab)\2p-data\2015-03-05_AF11\Resaved';
        session_nr=2;
        file_name=fullfile(data_folder,sprintf('20150305_AF11_%03d_resaved.tif',session_nr));
    case 2
        session_nr=1;
        data_folder='D:\Dropbox (coxlab)\2p-data\2015-08-03_AH02_init';
        file_name=fullfile(data_folder,sprintf('2015-08-03_AH02_%03d_resaved.tif',session_nr))   
end
%file_name='C:\Users\labuser\Desktop\20150305_AF11_009_resaved.tif';



info=imfinfo(file_name);

nFrames=length(info);
W=info(1).Width;
H=info(1).Height;
dataMatrix=zeros(nFrames,7);

offset=double(intmax('uint16')/2);
for iFrame=50%1:nFrames %min([500 nFrames])
    % read out last line of the frame, here we store frame specific
    % info: bitCodes, position, laser power
    frame=imread(file_name,iFrame,'info',info);
    imshow(frame,[])
    
    flyback_line=double(imread(file_name,iFrame,'info',info,'PixelRegion',{[H H],[1 W]}));
    parse_flyback_line(flyback_line)
end
%%
%plot(flyback_line)