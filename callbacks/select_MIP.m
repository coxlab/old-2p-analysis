function select_MIP(varargin)

H=varargin{1};
handles=guidata(H);
MIP_type=get(H,'value');

session_data=handles.session_data;
switch MIP_type
    case 1
        handles.MIP=session_data.MIP_avg.data;
        handles.MIP_gamma_val=session_data.MIP_avg.gamma_val;
    case 2
        handles.MIP=session_data.MIP_max.data;
        handles.MIP_gamma_val=session_data.MIP_max.gamma_val;
    case 3
        handles.MIP=session_data.MIP_std.data;
        handles.MIP_gamma_val=session_data.MIP_std.gamma_val;
    case 4
        handles.MIP=session_data.MIP_cc_local.data; 
        handles.MIP_gamma_val=session_data.MIP_cc_local.gamma_val;
end
handles.MIP_type=MIP_type;
handles.MIP=handles.MIP/max(handles.MIP(:))*256*256;
guidata(H,handles)

update_GUI(H)

