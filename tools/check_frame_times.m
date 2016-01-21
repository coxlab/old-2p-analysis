%clear all
%clc

%%
switch 1
    case 1
        file_name='C:\Users\labuser\Documents\ben\Movies\2015-09-28_dryRun_teamviewer\2015-09-28_dryRun004.tif';
end

info=imfinfo(file_name);

nFrames=length(info);
W=info(1).Width;
H=info(1).Height;
data_matrix=zeros(nFrames,7);

%%
offset=double(intmax('uint16')/2);
t0=clock;
for iFrame=1:nFrames %min([500 nFrames])
    % read out last line of the frame, here we store frame specific
    % info: bitCodes, position, laser power
    frame=imread(file_name,iFrame,'info',info);        
    flyback_line=double(imread(file_name,iFrame,'info',info,'PixelRegion',{[H H],[1 W]}));
    FBL=parse_flyback_line(flyback_line);
    data_matrix(iFrame,1:2)=[iFrame FBL.timestamp];
    progress(iFrame,nFrames,t0)
end

%%
T=data_matrix(:,2);
T=T-min(T);
IFI=diff(T);
frame_duration=mean(IFI);
frame_rate_est=1/frame_duration;
frame_rate=6.2;
frame_duration_actual=1/frame_rate;

figure(555)
plot(IFI)
hold on
plot([0 nFrames],[frame_duration_actual frame_duration_actual],'r-')
plot([0 nFrames],[frame_duration frame_duration],'g-')
hold off


%plot(flyback_line)