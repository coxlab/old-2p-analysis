function writeTable(varargin)

hObject=varargin{1};

if nargin>=2
    struct_name=varargin{2};
else
    error('Function requires at least 2 inputs...')
end
if nargin>=3
    index=varargin{3};
else
    index=1;
end

handles=guidata(hObject);
command=['S=handles.' struct_name '(' num2str(index) ');'];
eval(command)

parameters=fieldnames(S);
%eval(['parameters=fieldnames(handles.' struct_name ');']);
N=length(parameters);
for iParam=1:N    
    param_str=parameters{iParam};    
    param_str(param_str=='_')=' '; % convert to screen name, this will not affect the variable name used in handles
    
    command=['data{' num2str(iParam) ',1}=param_str;'];
    eval(command)    
    command=['data{' num2str(iParam) ',2}=S.' parameters{iParam} ';'];
    eval(command)
end

%data
set(hObject,'Data',data)


