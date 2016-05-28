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
        host=host(1:end-1);
        [~,user]=system('whoami');
        user=user(1:end-1);
    else
        %%% Linux
        platform='Unix';
        [~,host]=system('hostname');
        host=host(1:end-1);
        [~,user]=system('whoami');
        user=user(1:end-1);
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
