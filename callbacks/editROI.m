function editROI(varargin)
H=varargin{1};
handles=guidata(H);

% change alignment of the ROI without having to delete and create new,
% identity must be preserved
% show same area, but allow new points to be collected...

% select part from the MIP based on current ROI coords
ROI=handles.ROI(handles.ROI_selector);
ROI_rect=ROI.ROI_rect;

detail=handles.MIP(ROI_rect(2):ROI_rect(4),ROI_rect(1):ROI_rect(3));

FFT_mask=drawCircle2(size(detail,1)/4,[0 0],size(detail));
DC=mean(detail(:));
F=fft2(detail-DC);
detail=real(ifft2(fftshift(fftshift(F).*FFT_mask)))+DC;

%%% Prepare detail image to work on
handles.detail=detail;
handles.detail_gamma_val=1;

%%% Setup placeholder for new ROI
ROI_temp=handles.ROI_empty;
ROI_temp.base_coord=round(ROI.center_coords); % use actual center as new base_coord
ROI_temp.ROI_nr=ROI.ROI_nr;
handles.ROI_temp=ROI_temp;

handles.status=2;

%%% Show
set(handles.subplots(1).p(2),'xData',[],'yData',[])
set(handles.subplots(2).fig,'cLim',[min(handles.detail(:)) max(handles.detail(:))]);
set(handles.subplots(2).h(1),'cData',handles.detail);
set(handles.subplots(2).p(1),'xData',[],'yData',[])
set(handles.subplots(2).p(2),'xData',[],'yData',[])

guidata(H,handles)