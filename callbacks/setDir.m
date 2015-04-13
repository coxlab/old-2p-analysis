function setDir(varargin)
H=varargin{1};
handles=guidata(H);

data_folder=uigetdir(handles.data_root);

loadName=fullfile(data_folder,'data_analysis','session_overview.mat');

if exist(loadName,'file')==2
       
    load(loadName,'data_sessions')
    handles.nSessions=length(data_sessions);
        
    % populate session selector dropdown
    handles.SessionNames=strcat({'Session #'},num2str((1:handles.nSessions).')).';
    set(handles.session_selector_ui,'string',handles.SessionNames)
    set(handles.session_selector_ui,'value',1)
    handles.selected_session=get(handles.session_selector_ui,'value');
    
    handles.data_folder=data_folder;
    handles.data_sessions=data_sessions;
    
    guidata(H,handles)
    
    % load first session of folder
    loadData(H)
else
    disp('file not found')
    if isfield(handles,'data_folder')
        handles.data_folder
    end
end
