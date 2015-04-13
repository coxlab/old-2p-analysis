function keyUpFcn(varargin)
H=varargin{1};
key_event=varargin{2};
handles=guidata(H);

%key_event.Key
switch key_event.Key
    % modifiers
    case 'shift'
        handles.modifiers(1)=0;
    case 'control'
        handles.modifiers(2)=0;
    case 'alt'
        handles.modifiers(3)=0;
    case '0'  % command
        handles.modifiers(4)=0;
        
        
    case 'space'
        %handles.ROI
    otherwise
end


guidata(H,handles)