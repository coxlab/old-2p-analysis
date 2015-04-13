function selectROI(varargin)
H=varargin{1};
handles=guidata(H);

selected_ROI=1;
if nargin==2
    selected_ROI=varargin{2};
end

if nargin>=3
    selected_ROI=varargin{3};
end

handles.ROI_selector=selected_ROI;
handles.status=0;
guidata(H,handles)

update_GUI(H)
