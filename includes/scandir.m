function list = scandir(folder,varargin)

fullpath=0;
if nargin==1
    filetype='.*';
elseif nargin==2
    filetype=varargin{1};
elseif nargin==3
    filetype=varargin{1};
    fullpath=varargin{2};
end

cwd=cd;
cd(folder)
A=dir(['*' filetype]);
list=struct;
if isempty(A)
    list=[];
else
    for i=1:length(A)
        if fullpath==1
            list(i).name=fullfile(folder,A(i).name);
        else
            list(i).name=A(i).name;
        end
        list(i).date=A(i).date;
    end
end
cd(cwd)