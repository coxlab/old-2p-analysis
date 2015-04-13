function MIP_gamma(varargin)
H=varargin{1};
handles=guidata(H);
direction=varargin{3};

step_size=.05;
switch direction
    case 'up'
        handles.MIP_gamma_val=handles.MIP_gamma_val+step_size;
    case 'down'
        handles.MIP_gamma_val=handles.MIP_gamma_val-step_size;
end

%%% Keep track of this change in the datafile, so the gamma value can be
%%% reloaded
session_data=handles.session_data;
switch handles.MIP_type
    case 1        
        session_data.MIP_avg.gamma_val=handles.MIP_gamma_val;
    case 2
        session_data.MIP_max.gamma_val=handles.MIP_gamma_val;
    case 3
        session_data.MIP_std.gamma_val=handles.MIP_gamma_val;
    case 4
        session_data.MIP_cc_local.gamma_val=handles.MIP_gamma_val;
end
handles.session_data=session_data;

guidata(H,handles)

update_GUI(H)