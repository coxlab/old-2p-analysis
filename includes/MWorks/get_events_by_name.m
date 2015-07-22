function varargout=get_events_by_name(varargin)
%[events,nEvents,event_code]=get_events_by_name(mwk_file_name,tag_find [,event_codec])

%%% Parse inputs
mwk_file_name=varargin{1};
tag_find=varargin{2};

if nargin>=3&&~isempty(varargin{3})
    event_codec=varargin{3};
else
    A=getCodecs(mwk_file_name);
    event_codec=A.codec;
end

if nargin>=4&&~isempty(varargin{4})
    start_time=varargin{4};
else
    start_time=[];
end

if nargin>=5&&~isempty(varargin{5})
    end_time=varargin{5};
else
    end_time=[];
end

%%% build tag list and find code needed for selection
tag_list={event_codec.tagname}';
code_list=cat(1,event_codec.code);
found=ismember(tag_list,tag_find);
if any(found)
    code_selection=code_list(found);
    %%% Read events with given event_code from file
    if ~isempty(start_time)&&~isempty(end_time)
        events=getEvents(mwk_file_name, code_selection, start_time, end_time);
    else
        events=getEvents(mwk_file_name, code_selection);
    end
    
    %%% Construct outputs
    varargout{1}=events;
    varargout{2}=length(events);
    varargout{3}=code_selection;
else % No such tag
    %%% Construct outputs
    varargout{1}=[];
    varargout{2}=0;
    varargout{3}=-1;
end

