function varargout=get_events_by_name(varargin)

if ismac
    switch nargout
        case 1
            varargout{1}=get_events_by_name_mac(varargin{:});
        case 2
            [varargout{1},varargout{2}]=get_events_by_name_mac(varargin{:});
        case 3
            [varargout{1},varargout{2},varargout{3}]=get_events_by_name_mac(varargin{:});
    end
else
    switch nargout
        case 1
            varargout{1}=get_events_by_name_xplatform(varargin{:});
        case 2
            [varargout{1},varargout{2}]=get_events_by_name_xplatform(varargin{:});
        case 3
            [varargout{1},varargout{2},varargout{3}]=get_events_by_name_xplatform(varargin{:});
    end
end
