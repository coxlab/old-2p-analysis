function varargout=get_events_by_name_xplatform(varargin)
%[events,nEvents,event_code]=get_events_by_name_xplatform(mwk_file_name,tag_find [,event_codec])

%%% Parse inputs
mwk_file_name=varargin{1};
tag_find=varargin{2};

mat_file_name=strrep(mwk_file_name,'.mwk','.mat');
if exist(mat_file_name,'file')==2
    load(mat_file_name,'codecs','event_code_strings','event_code_selection','events')
else
    error('No converted .mwk file found for this session...')
end

if nargin>=3&&~isempty(varargin{3})
    %event_codec=varargin{3};
else
    %event_codec=codecs.codec;
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

%%% Get stim update events
[found,b]=ismember(tag_find,event_code_strings);
if found==1
    code_selection=cat(1,events.event_code)==event_code_selection(b);
    if ~isempty(start_time)&&~isempty(end_time)
        timestamp_vector=cat(1,events.time_us);
        code_selection=code_selection==1&between(timestamp_vector,[start_time end_time]);
    end
    MW_events=events(code_selection);
    
    %%% Construct outputs
    varargout{1}=MW_events;
    varargout{2}=length(MW_events);
    varargout{3}=code_selection;
else % No such tag
    %%% Construct outputs
    varargout{1}=[];
    varargout{2}=0;
    varargout{3}=-1;
end
