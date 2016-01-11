clear all
clc

header_script

expt = frGetExpt(exp_name);

folder_names{1}='151218_KS154_etl_2P_KS/run03_ori8_reversed_plane01';
folder_names{2}='151218_KS154_etl_2P_KS/run03_ori8_reversed2_plane02';
folder_names{3}='151218_KS154_etl_2P_KS/run03_movingbars_cardinal_plane01';
folder_names{4}='151218_KS154_etl_2P_KS/run03_movingbars_cardinal_plane02';
nFolders=length(folder_names);

%compare_vector=[1 2];
compare_vector=[3 4];

%% Load ROI definitions
load_name=fullfile(expt.dirs.analysis,folder_names{compare_vector(1)},'ROI_masks_ben.mat');
load(load_name)
ROI_props_01=ROI_props;
I1=CC_2D;
C1=coloredLabels;

load_name=fullfile(expt.dirs.analysis,folder_names{compare_vector(2)},'ROI_masks_ben.mat');
load(load_name)
ROI_props_02=ROI_props;
I2=CC_2D;
C2=coloredLabels;

%% check offsets
c=corr([I1(:) I2(:)]);
cc=normxcorr2(I1,I2);
subplot(211)
imshow(I1,[])
subplot(212)
imshow(I2,[])
colormap(green)

%% calc offsets
imshow(cc,[])
[r,c]=find(cc==max(cc(:)));
x=c-size(C1,2);
y=r-size(C1,1);

%% offset B
C2_offset=offsetIm(C2,-x,-y,0);
imshow(C2_offset,[])

c=corr(double([C1(:) C2_offset(:)]));
% [r,c]=find(cc==max(cc(:)));
% x=c-size(C1,2);
% y=r-size(C1,1);


%% look at data
subplot(311)
imshow(C1)
subplot(312)
imshow(C2)
subplot(313)
%imshow(mean(C1,3)+mean(C2,3)*2,[])
imshow(mean(C1,3)+mean(C2_offset,3)*2,[])


%% Compare individual ROIs for overlap
N1=length(ROI_props_01);
N2=length(ROI_props_02);

match_matrix=[];
match_index=1;
for iROI_01=1:N1
    M1=ROI_props_01(iROI_01).mask;
    for iROI_02=1:N2
        %[iROI_01 iROI_02]
        M2=ROI_props_02(iROI_02).mask;
        M2_offset=offsetIm(full(M2),-x,-y,0);
        if 0
            %%
            subplot(121)
            imshow(full(M1),[])
            subplot(122)
            imshow(M2_offset+full(M1),[])            
        end
        
        %% calc overlap
        %sum_im=full(M1)+full(M2);
        sum_im=full(M1)+full(M2_offset);
        %imshow(sum_im,[])
        overlap=sum(sum_im(:)==2)/sum(M1(:));
        
        if overlap>0
            %% register matching ROIs
            match_matrix(match_index,:)=[match_index iROI_01 iROI_02 overlap];            
            match_index=match_index+1;
        end        
    end
end

%%
match_matrix

for iMatch=1:size(match_matrix,1)
     A(:,:,iMatch)=full(ROI_props_01(match_matrix(iMatch,2)).mask);
     B(:,:,iMatch)=full(ROI_props_02(match_matrix(iMatch,3)).mask);
end
subplot(211)
imshow(mean(A,3),[])
subplot(212)
imshow(mean(B,3),[])


