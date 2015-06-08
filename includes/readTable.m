function readTable(varargin)

hObject=varargin{1};
Event=varargin{2};

if nargin>=3
    struct_name=varargin{3};
else
    error('Function requires at least 3 inputs...')
end
if nargin>=4
    index=varargin{4};
else
    index=1;
end

%%% Find row nr for cell that was changed
iParam=Event.Indices(1);

%%% Get handles and read field for requested properties
handles=guidata(hObject);
command=['S=handles.' struct_name '(' num2str(index) ');'];
eval(command)
parameters=fieldnames(S);

%%% Assign new value to selected parameter
new_value=Event.NewData;
if ischar(new_value) % handle chars
    command=['S.' parameters{iParam} '=''' new_value ''';'];
elseif islogical(new_value) % handle logics
    logical_strings={'false','true'};
    command=['S.' parameters{iParam} '=' logical_strings{new_value+1} ';'];
else % treat as double
    command=['S.' parameters{iParam} '=' num2str(new_value) ';'];
end
eval(command)


%%% Store S to handles
command=['handles.' struct_name '(' num2str(index) ')=S;'];
eval(command)

%%% Read all values and store parameter change
guidata(hObject,handles)