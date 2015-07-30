clear all
clc

load_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-07-20_AG02/2015-07-20_AG02_003.tif';

info=imfinfo(load_name); % never save info, is huge
H=info.Height;
W=info.Width;
cur_frame=double(imread(load_name,567,'info',info,'PixelRegion',{[1 H-1],[1 W]}));

ref=offsetIm(cur_frame,30,40);

%%% Calc normxcorr2
[CC_max,offset]=im_align(cur_frame,ref);

%%% Calc phase shift
[r,c]=PCdemo(cur_frame,ref);

%%

M=session_data.motion_correction.shift_matrix;

subplot(211)
plot(M(:,2))
hold on
plot(M(:,4))
hold off

subplot(212)
plot(M(:,3))
hold on
plot(M(:,5))
hold off

%%
session_data.visualize_motion_correction()