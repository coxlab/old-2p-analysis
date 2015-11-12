function shiftROI(varargin)
H=varargin{1};
mode=varargin{3};
direction=varargin{4};
handles=guidata(H);

step_size=1;
switch direction
    case 'up'
        offset=[0 step_size];
    case 'left'
        offset=[step_size 0];
    case 'right'
        offset=[-step_size 0];
    case 'down'
        offset=[0 -step_size];
end

% move all coordinates with the offset found
% 2DO: check if coords are valid after move!
switch mode
    case 'all'
        ROI=handles.ROI;
        nROI=length(ROI);
    case 'single'
        ROI=handles.ROI(handles.ROI_selector);
        nROI=length(ROI);        
end

for iROI=1:nROI
    ROI(iROI).base_coord=ROI(iROI).base_coord-offset;
    ROI(iROI).coords=ROI(iROI).coords-repmat(offset,size(ROI(iROI).coords,1),1);
    ROI(iROI).ellipse_coords=ROI(iROI).ellipse_coords-repmat(offset,size(ROI(iROI).ellipse_coords,1),1);
    ROI(iROI).coords_MIP=ROI(iROI).coords_MIP-repmat(offset,size(ROI(iROI).coords_MIP,1),1);
    %ROI(iROI).coords_MIP_plot=ROI(iROI).coords_MIP_plot-repmat(offset,size(ROI(iROI).coords_MIP_plot,1),1);
    ROI(iROI).center_coords=ROI(iROI).center_coords-repmat(offset,size(ROI(iROI).center_coords,1),1);
    ROI(iROI).ellipse_coords_centered=ROI(iROI).ellipse_coords_centered-repmat(offset,size(ROI(iROI).ellipse_coords_centered,1),1);
    ROI(iROI).ROI_rect=ROI(iROI).ROI_rect-repmat(offset,size(ROI(iROI).ROI_rect,1),2);
end
%fprintf('Shifted all coordinates by x=%d and y=%d \n',offset)


%%% Check if coords are still valid after the shift
if isfield(handles.session_data,'data')
    Height=handles.session_data.data(4);
    Width=handles.session_data.data(3);
else
    Height=handles.session_data.mov_info.Height;
    Width=handles.session_data.mov_info.Width;
end

if 0
    remove_vector=zeros(nROI,1);
    for iROI=1:nROI
        ROI_rect=ROI(iROI).ROI_rect;
        if any(~between(ROI_rect([1 3]),[1 Width]))||any(~between(ROI_rect([2 4]),[1 Height]))
            %disp('ROI is invalid after shift, removing')
            remove_vector(iROI)=1;
        else
            %disp('ROI still valid')
        end
    end
    ROI(remove_vector==1)=[];
    nRemoved=sum(remove_vector);
    if nRemoved>0
        fprintf('Removed %d ROI that became invalid after shift \n',nRemoved)
        handles.ROI_selector=1;
    end
end
%%% Copy edited ROI struct back to main data
switch mode
    case 'all'
        handles.ROI=ROI;
        %handles.ROI_selector=1;
        handles.nROI=length(handles.ROI);
        handles.ROI_temp=handles.ROI_empty;
    case 'single'
        %if nRemoved==0
            handles.ROI(handles.ROI_selector)=ROI;
        %else
%             handles.ROI(handles.ROI_selector)=[];
%             handles.nROI=length(handles.ROI);
%             handles.ROI_temp=handles.ROI_empty;
%             disp('we lost the ROI')
%         end
end

guidata(H,handles)

update_GUI(H)