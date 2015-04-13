function selectSession(varargin)
H=varargin{1};
handles=guidata(H);

handles.selected_session=get(H,'value');

guidata(H,handles)

loadData(H)