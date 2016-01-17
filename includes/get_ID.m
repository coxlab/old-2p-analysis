function varargout=get_ID(varargin)
%%% Returns: hostname, username and platform
% Cross-platform function to get host and user names

if ispc
    %%% Windows
    platform='Windows';
    host=getenv('computername');
    user=getenv('username');
elseif isunix
    if ismac
        platform='Mac';
        %%% Mac
        [~,host]=system('hostname');
        [~,user]=system('whoami');
    else
        %%% Linux
        platform='Unix';
        [~,host]=system('hostname');
        [~,user]=system('whoami');
    end
else
    %% Unknown
    platform=[];
    host=[];
    user=[];
end

varargout{1}=host;
varargout{2}=user;
varargout{3}=platform;
