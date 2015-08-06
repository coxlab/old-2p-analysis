function varargout=parse_flyback_line(varargin)

if nargin>=1
    vector=varargin{1};
else
    error('No input')
end

%%% Find all bitCodes
if nargin>=2&&~isempty(varargin{2})
    nBitCodes=varargin{2};
else
    nBitCodes=find(vector==0,1,'first')-1;
    if nBitCodes==0
        nBitCodes=[];
    end
end

%% Determine version of flyback format
if vector(1)>20 % newest
    version=vector(1);
else
    if vector(end-17)>2000 % current
        version=10;
    else % old
        version=1;
    end
end

offset=double(intmax('uint16')/2);

switch version
    case 21
        nBitCodes=vector(2);
        
        %%% handle bitcodes
        bitCode_vector=vector(1:nBitCodes);
        main_bitCode=mode(bitCode_vector);
        
    case 10
        %%% handle bitcodes
        bitCode_vector=vector(1:nBitCodes);
        main_bitCode=mode(bitCode_vector);
        
        %%% handle headers
        headers=vector(end,end-17:end);
        date_num=headers(1:6)./[1 1e3 1e3 1e3 1e3 1e3]; % time in days
        timestamp=datenum(date_num)*24*60*60; % in seconds
        switch_times=headers(7:10); % ms
        switch_detected=any(switch_times>0);
        xyz_hi_res=(headers(11:13)-offset)/10; % 0.1 micron
        xyz=headers(14:16)-offset; % micron
        piezo=headers(17); % micron
        laser_power=headers(18); % percent! not mW
        
        results=struct;
        results.bitCode_vector=bitCode_vector(:);
        results.main_bitCode=main_bitCode;
        results.nBitCodes=nBitCodes;
        results.date_num=date_num;
        results.timestamp=timestamp;
        results.switch_times=switch_times;
        results.switch_detected=switch_detected;
        results.xyz_micron=xyz;
        results.xyz_submicron=xyz_hi_res;
        results.piezo=piezo;
        results.laser_power=laser_power;
    case 1
        %%% handle bitcodes
        bitCode_vector=vector(1:nBitCodes);
        main_bitCode=mode(bitCode_vector);
        
        headers=vector(end,end-4:end);
        
        xyz=headers(1:3)-offset;
        piezo=headers(4); % micron
        laser_power=headers(5); % percent! not mW
        
        results=struct;
        results.bitCode_vector=bitCode_vector(:);
        results.main_bitCode=main_bitCode;
        results.nBitCodes=nBitCodes;
        
        results.date_num=[];
        results.timestamp=[];
        results.switch_times=[];
        results.switch_detected=[];
        results.xyz_micron=xyz;
        results.xyz_submicron=[];
                
        results.piezo=piezo;
        results.laser_power=laser_power;
end

varargout{1}=results;
varargout{2}=nBitCodes;


