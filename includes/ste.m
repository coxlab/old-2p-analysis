function varargout = ste(varargin)

vector=varargin{1};
%if any(size(vector)==1)
%    vector=vector(:);
%end
if nargin>=2
    dim=varargin{2};
else
    dim=1;
end
N=length(vector);
stdev=std(vector,[],dim);

varargout{1}=stdev/sqrt(N);