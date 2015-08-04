clear all
clc

file_name='C:\Users\labuser\Desktop\2015-08-03_AH02_012.tif';

info=imfinfo(file_name);

nFrames=length(info);
W=info(1).Width;
H=info(1).Height;



dataMatrix=zeros(nFrames,7);

offset=double(intmax('uint16')/2);
for iFrame=1:nFrames %min([500 nFrames])
    % read out last line of the frame, here we store frame specific
    % info: bitCodes, position, laser power
    flyback_line=double(imread(file_name,iFrame,'info',info,'PixelRegion',{[H H],[1 W]}))
end
%%
plot(flyback_line)