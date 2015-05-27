function varargout=parse_conditions(varargin)

X=varargin{1};

diff([-1;X;-1])

varargout{1}=0;