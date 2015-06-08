function writeSummaryTable(varargin)

hObject=varargin{1};

if nargin>=2
    struct_name=varargin{2};
else
    error('Function requires at least 2 inputs...')
end
% if nargin>=3
%     index=varargin{3};
% else
%     index=1;
% end

handles=guidata(hObject);
%command=['S=handles.' struct_name '(' num2str(index) ');'];
command=['S=handles.' struct_name ';'];
eval(command)

parameters=fieldnames(S);
%eval(['parameters=fieldnames(handles.' struct_name ');']);
nRecords=length(S);
nParams=length(parameters);
for iParam=1:nParams
    param_str=parameters{iParam};
    param_str(param_str=='_')=' '; % convert to screen name, this will not affect the variable name used in handles
    Columnnames{iParam}=param_str;
    
    for iRecord=1:nRecords
        
        %command=['data{' num2str(iParam) ',1}=param_str;'];
        %eval(command)
        command=['data{' num2str(iRecord) ',' num2str(iParam) '}=S(' num2str(iRecord) ').' parameters{iParam} ';'];
        eval(command)
    end
end

set(hObject,'Columnname',Columnnames)
%data
set(hObject,'Data',data)


