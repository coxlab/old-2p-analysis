function saveData(varargin)
H=varargin{1};
handles=guidata(H);
session_data=handles.session_data;
if isfield(session_data.ROI_definitions,'ROI_nr')
    %session_data.ROI_definitions
    % Clear old ROI_definitions struct
    session_data=rmfield(session_data,'ROI_definitions');
    disp('Cleared old ROI_definitions struct...')
end

session_data.ROI_definitions(handles.ROI_definition_nr).ROI=handles.ROI;
saveName=handles.loadName;
save(saveName,'session_data')
disp(['Data saved to ' saveName])