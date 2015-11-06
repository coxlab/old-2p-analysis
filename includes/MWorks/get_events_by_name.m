function out=get_events_by_name(varargin)

if ismac
    out=get_events_by_name_mac(varargin{:});
else
    out=get_events_by_name_xplatform(varargin{:});
end