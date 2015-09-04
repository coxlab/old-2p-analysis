clear all
%clc

%file_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-05_AF11/20150305_AF11_009_retagged.tif';
file_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-05_AF11/20150305_AF11_009_resaved.tif';
%file_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-03_AH02_init/2015-08-03_AH02_001.tif';


obj=Tiff(file_name,'r');
%obj.getTagNames()
%session_data.import_movie(file_name)
%%
%obj.getTag('ImageDescription')
%obj.setTag('ImageDescription','Removed tag...')
%obj.getTag('ImageDescription')

%%
im=double(obj.read);

imagesc(corr(im'))

%%
obj.TagID


%% Go deeper using tifflib

FileID = tifflib('open',file_name,'r');

InfoImage=imfinfo(file_name);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);

tifflib('setDirectory',FileID,0);
rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);

for r = 1:rps:nImage
    row_inds = r:min(nImage,r+rps-1);
    stripNum = tifflib('computeStrip',FileID,r);
    
    S=tifflib('readEncodedStrip',FileID,stripNum-1);
end

%%
tifflib('setDirectory',FileID,2);
S1=tifflib('readEncodedStrip',FileID,5);
S2=tifflib('readEncodedStrip',FileID,6);

M=double([S1' S2']);
imagesc(corr(M))
axis square


%%
tifflib('retrieveMetadata',FileID)
