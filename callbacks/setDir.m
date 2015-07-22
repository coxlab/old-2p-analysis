function setDir(varargin)
H=varargin{1};
handles=guidata(H);

data_folder=uigetdir(handles.data_root);

%loadName=fullfile(data_folder,'data_analysis','session_overview.mat');
mat_folder=fullfile(data_folder,'data_analysis');

%if exist(loadName,'file')==2
if isdir(mat_folder)
    %load(loadName,'data_sessions')
    data_sessions=scandir(mat_folder,'mat');
    handles.nSessions=length(data_sessions);
        
    % populate session selector dropdown
    %handles.SessionNames=strcat({'Session #'},num2str((1:handles.nSessions).')).';
    handles.SessionNames={data_sessions.name};
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
