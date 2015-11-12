ROIs=get_ROI_definitions(session_data);
N=length(ROIs);

ROI_matrix=struct;

for iROI=1:N
    FOV_rect_px=session_data.FOV_info.size_px;
    FOV_rect_um=session_data.FOV_info.size_um;
    
    ROI=ROIs(iROI);
    
    %%% Construct mask for soma
    mask_soma=poly2mask(ROI.coords_MIP(:,1),ROI.coords_MIP(:,2),session_data.FOV_info.size_px(2),session_data.FOV_info.size_px(1));
    
    stretch_factor=round(session_data.FOV_info.pixel_size_micron/min(session_data.FOV_info.pixel_size_micron));
    stretch_coords=round(fliplr(session_data.FOV_info.size_px).*stretch_factor);
    
    %%% Construct mask for neuropil    
    mask_soma=imresize(mask_soma,stretch_coords,'nearest');
    %mask_neuropil=imdilate(mask_soma,se)-mask_soma;
    mask_neuropil=bwmorph(mask_soma,'thicken',4)-mask_soma;
    
    mask_soma=imresize(mask_soma,fliplr(FOV_rect_px),'nearest');
    mask_neuropil=imresize(mask_neuropil,fliplr(FOV_rect_px),'nearest');
    
    ROI_matrix(iROI).mask_soma=sparse(mask_soma);
    ROI_matrix(iROI).mask_neuropil=sparse(mask_neuropil);
    ROI_matrix(iROI).data=[iROI ROI.ROI_nr sum(mask_neuropil(:)) sum(mask_soma(:)) sum(mask_neuropil(:))/sum(mask_soma(:))];
    
    P=regionprops(mask_soma,'boundingBox');
    ROI_matrix(iROI).rect_soma=round(P.BoundingBox);
    P=regionprops(mask_neuropil,'boundingBox');
    ROI_matrix(iROI).rect_neuropil=round(P.BoundingBox);
end

%%
data=cat(1,ROI_matrix.data);
mean(data(:,end))
mask_soma=combine_sparse_masks(ROI_matrix,'mask_soma');
mask_neuropil=combine_sparse_masks(ROI_matrix,'mask_neuropil');

%%
%mask=sum(ROI_matrix_full,3);

figure(132)
imshow(imresize(mask_neuropil*.5+mask_soma,stretch_coords,'nearest'),[])

%% find overlap between neuropil regions and soma regions
iROI=N;
hole=ROI_matrix(iROI).mask_neuropil+mask_soma>1;
A=ROI_matrix(iROI).mask_neuropil;
B=A;
B(hole==1)=0;
imshow(full(A)*.5+full(B),[])

