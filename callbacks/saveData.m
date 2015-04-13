function saveData(varargin)
H=varargin{1};
handles=guidata(H);
session_data=handles.session_data;
session_data.ROI_definitions=handles.ROI;
saveName=handles.loadName;
save(saveName,'session_data')
disp(['Data saved to ' saveName])