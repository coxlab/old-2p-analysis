clear all
clc

%2DO:
% V add FOV size and position
% V do ROI stuff => other script
% - run script on all data you have MWorks events for!
% V linking SI and MW has been finalized, check upon generalization though!

% found error in bitCode matching for session 2 of /Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-05_AF03 row 178 of bitCode matched matrix
% turns out that acquisition was hanging and many bitCodes were stored into
% the flyback line... hard to fix... the main samplesAcquired function will
% not have saved frames, while the miniloop checking the pixelclock did
% accumulate, so this captured the actual frames presented, but we have no
% corresponding 2p data for this...
% sort of got around it
% found another error in 2015-03-04_AF11, session 4: bitcodes are off as of
% frame 201, but come together some frames later. Cut those out for now.

% consider saving frame time as clock into six bytes : uint16([2015 12 31 24 60 60].*[1 ones(1,5)+1000])

%   V check feasibility of other frame matching method based on difference between MW events, estimate number of frames presented.
%       V check using actually properly processed file, all good files
%       match up perfectly!
%       => using it on bogus flyback line session

% check super resolution from debug files

% BV20150407: Found session 02 of 04007_AF11 to be bogus, bitcodes are only in first
% pixel of the flyback line frames 1342-1344

% BV20150408: need a way to distinguish between experiment types. This
% should be extractable from the mworks file, <! find out how !>

% BV20150420: need to make -nojvm version where we supply a folder to the
% script as a function so we can run on the server, if not on linux use uigetdir

header_script

%%
%%% Use uigetdir for general use
%cd(data_root)
%data_folder=uigetdir(data_root);

%%% Gather all movies
files=scandir(data_folder,'tif');
nFiles=length(files);

%%% Prepare MWorks analysis part
mwk_files=scandir(data_folder,'mwk');
if length(mwk_files)>1
    error('Found more MWorks datafiles then expected... Make sure each MWorks file is located in a separate folder together with the corresponding scanimage movies')
elseif isempty(mwk_files)
    error('No valid .mwk file found')
else
    mwk_file_name=fullfile(data_folder,mwk_files(1).name); % should be 1 file per folder
end

% %%% Check existence and type of .mwk file
if exist(mwk_file_name,'file')==2 % found file
elseif exist(mwk_file_name,'file')==7 % found dir
    % file can be read whehter it is in it's folder of in main folder
else
    error('No valid .mwk file found')
end

%%
fprintf(['Processing %d files from folder' data_folder '\n'],nFiles)
t0=clock;
file_info=struct;
field_names={'Index','nFrames','Width','Height','FrameRate','StartTime','std(X)','std(Y)','std(Z)','DataType'};
for iFile=1:nFiles
    file_name=fullfile(data_folder,files(iFile).name);
    info=imfinfo(file_name);
    
    %%% Get basic info about movie
    nFrames=length(info);
    W=info(1).Width;
    H=info(1).Height;
    
    %% Get SCIM headers
    scim_info=info(1).ImageDescription;
    scinfo=strsplit(char(13),scim_info);
    scim_info=strjoin(scinfo,[';' char(13)]);
    eval([scim_info  ';']); % revive scan image variables
    frame_rate=state.acq.frameRate;
    frame_start_time=state.internal.softTriggerTimeString;
    
    %%% Extract position, depth and laser power information
    % distinguish between stacks and ROI data_types
    dataMatrix=zeros(nFrames,7);
    
    offset=double(intmax('uint16')/2);
    for iFrame=1:nFrames %min([500 nFrames])
        % read out last line of the frame, here we store frame specific
        % info: bitCodes, position, laser power
        flyback_line=double(imread(file_name,iFrame,'info',info,'PixelRegion',{[H H],[1 W]}));
        
        
        %%% read and interpret bitCodes, still crude
        switch H
            case 512
                bitCode_length=64;
            case 400
                bitCode_length=50; % # of bitcode samples per acquis. frame
            case 200
                bitCode_length=0;
        end
        if strcmpi(data_folder,'/Users/benvermaercke/Dropbox (coxlab)_offline/2p-data/2015-02-27_debug')
            if iFile==3
                bitCode_length=0;
            end
        end
        
        bitCode=mode(flyback_line(1:bitCode_length)); % just pick the most frequently occuring bitcode in the sequence for that frame
        
        %%% Correcting error in second session of 0407_AF11
        if strcmpi(data_folder,'/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-07_AF11')
            if iFile==2
                if ismember(iFrame,1342:1344)
                    bitCode=flyback_line(1);
                end
            end
        end
        %flyback_line(1:bitCode_length+1)
        
        
        % There should always be 6 OR 7 frames for the "blanks" and always
        % 6 frames for the stimuli (depending on the frame rate and stim
        % duration:  2 sec on/off ~ 6 frames on/off @ 3 Hz)
        
        % flyback_line
        % should use the crossings in the future, not sure about correctness
        % yet: crossings=flyback_line(end-8:end-5);
        
        % Crossings keep track of the exact sample # relative to each sub-cycle
        % of samples, which can be used to reconstrut exact time of
        % threshold crossing (i.e., relative to sub-cycle, not the absolute
        % time of the frame start).
        
        %%% read position info and rescale if nec
        position_laser=flyback_line(end-4:end);
        if W==256
            % don't do it for the bogus files
        else
            if position_laser(1)-offset>-10000
                position_laser=position_laser-[offset offset offset 0 0];
            end
        end
        
        %%% build data matrix
        dataMatrix(iFrame,1:7)=[iFrame bitCode position_laser];
    end
    
    %%
    %%% Use these data to specify data type
    % options:
    % 0) undefined
    % 1) random, compromised flyback line
    % 2) scanning around, flyby
    % 3) FOV with small correction
    % 4) stack, trajectory GUI
    % 5) stack, manual
    % 6) FOV => match to MWorks
    data_type=0;
    X=dataMatrix(:,3);Y=dataMatrix(:,4);Z=dataMatrix(:,5);
    
    if strcmpi(data_folder,'/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-16_AF11')
        % messed with the coords during these sessions
        if iFile==3
            X=dataMatrix(:,3)*0+dataMatrix(1,3);
            Y=dataMatrix(:,4)*0+dataMatrix(1,4);
            Z=dataMatrix(:,5)*0+dataMatrix(1,5);
        end
        if iFile==4
            X=dataMatrix(:,3)*0+dataMatrix(1,3);
            Y=dataMatrix(:,4)*0+dataMatrix(1,4);
            sel=dataMatrix(:,3)==dataMatrix(1,3);
            Z=dataMatrix(:,5);
            Z(sel)=Z(sel)-49;
        end
    end
    
    TH_scanning=1000; % if std for x, y, and z coordinate vectors b/w 100-1000, case 2
    TH_stack=80; % if std less than this for x, y, and z, NOT case 4 (std should be greater than 100 for stacks)...
    TH_small_shift=50; % std for x, y, and z (case 3)
    TH_no_shift=5; % std for x, y, and z (case 6)
    
    % If changes in X or Y, bogus data, just scanning or small correction,
    % compare std's
    if std(X)>TH_scanning&&std(Y)>TH_scanning&&std(Z)>TH_scanning % probably no valid flyback line
        data_type=1;
    end
    if any(between([std(X) std(Y) std(Z)],[TH_small_shift TH_scanning]))
        %if std(X)<TH_scanning&&std(Y)<TH_scanning&&std(Z)<TH_scanning
        data_type=2;
    end
    
    if any(between([std(X) std(Y) std(Z)],[TH_no_shift TH_small_shift]))
        %if std(X)>TH_small_shift||std(Y)>TH_small_shift||std(Z)>TH_small_shift
        data_type=3;
    end
    
    % If changes in Z, doing stack, make sure changing part has fixed
    % interval
    TH_manual_stack=.5; % if std is higher than this, then it's a manul stack since w/ trajectoryGui speed is constant down z-axis
    if std(X)<TH_no_shift&&std(Y)<TH_no_shift&&std(Z)>TH_stack
        
        % detect changing part of Z coords
        % if changes are uniform 4 else 5
        change_vector.diff=between(diff(Z),[1 3]);
        change_vector.start=find(change_vector.diff==1,1,'first');
        change_vector.end=find(change_vector.diff==1,1,'last');
        if std(change_vector.diff(change_vector.start:change_vector.end))<TH_manual_stack
            data_type=4;
        else
            data_type=5;
        end
    end
    
    % No changes, FOV ROI exp
    if std(X)<TH_no_shift&&std(Y)<TH_no_shift&&std(Z)<TH_no_shift
        data_type=6;
    end
    
    if data_type==0
        error('No data type could be determined...')
    end
    
    %%% These numbers allow to get position and FOV size, so we can plot
    %%% these in window coordinates
    if data_type==6
        FOV_info.center=[mean(X) mean(Y)]; % micron
        FOV_info.z_depth=mean(Z); % micron
        FOV_info.size_px=[H W]; % pixels
        FOV_info.pixel_size_micron=.85;
        FOV_info.size_um=FOV_info.size_px*FOV_info.pixel_size_micron; % micron
    else
        FOV_info=struct;
    end
    
    %%% Apply medfilter to filter out glitches
    % be careful that this does not induce other strange behaviors!
    dataMatrix(:,2)=medfilt1(dataMatrix(:,2),3);
    
    %%% collect data per file
    file_info(iFile).file_name=file_name;
    file_info(iFile).expType=[];
    file_info(iFile).expType_name=[];    
    %file_info(iFile).info=info; % do not save, this is huge and can be
    % rapidly read out from the datafile at any point
    file_info(iFile).state=state;
    file_info(iFile).field_names=field_names;
    file_info(iFile).data=[iFile nFrames W H frame_rate datenum(frame_start_time) std(X) std(Y) std(Z) data_type];
    file_info(iFile).dataMatrix=dataMatrix;
    file_info(iFile).FOV_info=FOV_info;
    file_info(iFile).ROI_definitions=struct;
    
    progress(iFile,nFiles,t0)
end

if 0
    %%
    T=now;
    datestr(T)
    datestr(T+60/1000/60)
    
    %%
    T1='06-Mar-2015 15:55:00'
    T2='06-Mar-2015 15:56:00'
    minute=datenum(T2)-datenum(T1);
    minute*60*24 % 1440 conversion factor from now() units
end

%%
SI_session_data=cat(1,file_info.data);
SI_session_data(:,end)

%%% Select only session that have stable imaging of single FOV
sel=ismember(SI_session_data(:,10),[1 6]);
nExpSessions=sum(sel); % These need matching to MWorks

%%% Select start and end times for all selected session
SI_session_data=SI_session_data(sel,:);

SI_session_start_times=SI_session_data(:,6);
SI_session_start_times=SI_session_start_times-SI_session_start_times(1);
SI_session_start_times=SI_session_start_times*60*24; % unit of now() to minutes

SI_session_duration=SI_session_data(:,2)./SI_session_data(:,5)/60; % seconds to minutes

SI_session_times=[SI_session_start_times [diff(SI_session_start_times);SI_session_duration(end)]]; % no length value for last movie, since it has noone following it. Use session duration
SI_session_allocation=[SI_session_data(:,1) SI_session_times SI_session_duration SI_session_times(:,2)-SI_session_duration];


TH_session_duration=.7; % minute
SI_session_allocation(SI_session_allocation(:,4)<TH_session_duration,:)=[];








%% Read MWorks file and extract events
tic
disp('Reading MWK file...')
A=getCodecs(mwk_file_name);
event_codec=A.codec;

code_selection=8; % collect stim_display_update events
MW_events=getEvents(mwk_file_name, code_selection);
nEvents=length(MW_events);

fprintf('Loading MWorks events took %3.2f seconds.\n',toc)

if 0 % testing codes
    %%
    for iCodec=1:length(event_codec)
        disp([num2str(event_codec(iCodec).code) ': ' event_codec(iCodec).tagname])
    end
    
    code_selection=66; % collect stim_display_update events
    MW_events=getEvents(mwk_file_name, code_selection);
    nEvents=length(MW_events);
    
    %%
    for iEvent=1:nEvents
        disp(double([(MW_events(iEvent).time_us-MW_events(1).time_us)/1e6 MW_events(iEvent).data]))
    end
    
end

%% BV20150416: Define experiment type
% Right now, we have 2 possible experimental protocols, but this will
% increase with time. We need a unambiguous methods to label each. This
% expType variable will determine how MWorks events are processed.
% Perhaps we could add a variable expType to each protocol and read from
% the codec which field it got assigned to (this may vary depending on
% number of variables). Then read it out.
% Since we cannot change MWorks protocol and keep streaming to same file,
% assume 1 expType per file. Different protocols should then have different
% folders!

% For all old files, see if we have an expType field in the codec
tag_names={event_codec.tagname}';
code_names=cat(1,event_codec.code);
exp_tag='ExpType';
%exp_tag='stop_it';
tag_nr=find(ismember(tag_names,exp_tag),1);
expType=0;
%expType_name='';
expType_name_vector={'RSVP','Retinomapping'};
if ~isempty(tag_nr) % if expType tag is present, just read it out
    code_selection=code_names(tag_nr);
    MW_events_expType=getEvents(mwk_file_name, code_selection);
    nEvents=length(MW_events_expType);
    type_vector=cat(1,MW_events_expType.data);
    expType=mode(type_vector);
else % for older files, find work-around
    %%% Find tags specific to a certain experiment =risky
    tag_nr=find(ismember(tag_names,'stm_pos_x'),1);
    if ~isempty(tag_nr)
        expType=1; % RSVP
        %expType_name='RSVP';
    end
    
    %%% Find tags specific to a certain experiment =risky
    tag_nr=find(ismember(tag_names,'show_vertical_bar'),1);
    if ~isempty(tag_nr)
        expType=2; % Retinomapping
        %expType_name='Retinomapping';
    end
    
    if expType==0
        disp('Unable to determine experiment type, how do you wish to proceed?')
    end
end
expType_name=expType_name_vector{expType};

for iFile=1:nFiles
    file_info(iFile).expType=expType;
    file_info(iFile).expType_name=expType_name;
end

%% Set experiment specific variables
switch expType
    case 0
        error('No experiment type defined...')
    case 1
        stim_duration=[2.5 2]*1e6;
    case 2
        stim_duration=[4 .0167]*1e6;
end

%% Get timestamps
timestamps=double(cat(1,MW_events.time_us));
timestamps=timestamps-timestamps(1);
% Build timestamp matrix with starttime and duration of events
TS_matrix=[timestamps [diff(timestamps);0]];

% Find all breaks, longer than stim duration seconds
%max_stim_time=max([2.1 1])*1e6;
max_stim_time=max(stim_duration);
break_vector=TS_matrix(:,2)>max_stim_time;
nBreaks=sum(break_vector)

% These times indicate the idle time between session
start_break=find(break_vector);
end_break=start_break+1;

% calc number of event in session
nEvents_in_session=diff([[0;end_break] [start_break;0]],[],2)+1;

A=[TS_matrix(1,1);TS_matrix(end_break)]/1e6/60; % start of gap
B=[TS_matrix(start_break);TS_matrix(end,1)]/1e6/60; % end of gap

MW_session_times=[A B [diff(A);0] diff([A B],[],2) [diff(A);0]-diff([A B],[],2) A*60e6 B*60e6 [[0;end_break] [start_break;0]] nEvents_in_session];

%%% impose minimal duration for valid sessions
MW_session_times(MW_session_times(:,4)<TH_session_duration,:)=[];

%%% adhoc solution to get debug session processed.
if strcmpi(data_folder,'/Users/benvermaercke/Dropbox (coxlab)_offline/2p-data/2015-02-27_debug')
    MW_session_times(1,:)=[];
    file_info(3).data(end)=1; % create artificial bitCodes to replace faulty SI_bitCodes
end


% Rescale time to this first element of the cropped matrix
MW_session_times(:,2)=MW_session_times(:,1)-MW_session_times(1,1);

% Renumber
MW_session_times(:,1)=1:size(MW_session_times,1);


%% Try to match each recorded ROI file to a MWorks block of events
try
    % try to match automatically by looking at start of each session / trial
    % block
    A=SI_session_allocation(:,2);
    B=MW_session_times(:,2);
    
    TH_time_offset=.010; % in seconds
    match_vector=zeros(length(A),1);
    for iSess=1:length(A)
        [dist,pos]=min(abs(A(iSess)-B));
        if dist<TH_time_offset
            match_vector(iSess)=pos;
        else
            error('Could not find a match within the specified time interval, try manual matching')
        end
    end
    
    MW_session_times_matched=MW_session_times(match_vector,:);
    
    % compare session durations as a sanity check
    TH_duration_mismatch=0.500; % in seconds
    if any(abs(diff([SI_session_allocation(:,4) MW_session_times_matched(:,4)],[],2))>TH_duration_mismatch)
        error('Difference between session duration detected, go manual')
    end
    disp('Passed sanity check: session duration')
    
    matching_matrix=[SI_session_allocation(:,1) MW_session_times_matched(:,1)]
    
catch % failed, do it manually
    % use read and write tables we designed for the ephys data analysis
    % pipeline
    if strcmpi(data_folder,'/Users/benvermaercke/Dropbox (coxlab)_offline/2p-data/20150225_jrat3_ben')
        matching_matrix=[1 1 ; 2 2];
        MW_session_times_matched=MW_session_times;
        file_info(1).data(end)=1;
        file_info(2).data(end)=1;
    end
    
    E=lasterror;
    E.message
end


%%
%%% Assign stimulus properties to each frame of each movie
% - depending on whether bitCodes from SCIM files are useable use those or
% allocate based on times in between MWorks display update events
% 2DO everify whether this approach is valid in the case where we have all
% information

% we need for each session:
% SI:
% - number of frames
% - bitCodes

% MW:
% - start and stop timestamps
% - stim info
% - bitCodes

nSessions=size(matching_matrix,1);
data_sessions=struct('file_name',[],'expType',[],'expType_name',[],'state',struct,'field_names',[],'data',[],'dataMatrix',[],'FOV_info',struct,'ROI_definitions',struct,'stimulus_matrix_ext',[]);
good_sessions=0;
for iSess=1:nSessions
    F=file_info(matching_matrix(iSess,1));
    data_type=F.data(end);
    
    %% SI
    nFrames=F.data(:,2);
    frame_rate=F.data(5);
    bitCodes_scanimage=F.dataMatrix(:,1:2);
    
    switch expType
        case 1
            %% MW: collect stimulus information
            sess_events=MW_events(MW_session_times_matched(iSess,8):MW_session_times_matched(iSess,9));
            T=double(cat(1,sess_events.time_us));
            % first stim was triggered by the end of the first SI frame so
            % appears at second frame. first frame is always free of stimulus, plus
            % it has an artifact of the delayed shutter opening => ignore
            
            nSess_events=length(sess_events);
            count=1;
            stim_counter=1;
            stimulus_matrix=[];
            for iEvent=1:nSess_events
                D=sess_events(iEvent).data; % # of session events
                timestamp=double(sess_events(iEvent).time_us);
                
                if length(D)==3 % stim
                    D{2}.stim_nr=str2double(D{2}.name);
                    stimulus_matrix(count,:)=[timestamp 0 double(D{3}.bit_code) stim_counter D{2}.stim_nr D{2}.pos_x D{2}.pos_y];
                    stim_counter=stim_counter+1;
                    count=count+1;
                elseif length(D)==2 % blank
                    stimulus_matrix(count,:)=[timestamp 0 double(D{2}.bit_code) -1 -1 -1 -1];
                    count=count+1;
                else
                    % empty event, can occurs at begin and end of MW session file
                end
            end
            
            % add stimulus duration
            stimulus_matrix(:,2)=[diff(stimulus_matrix(:,1));0];
            
            %% link SI and MW events
            if data_type==6 % default FOV
                try
                    cur_bitCode=-1;
                    stim_properties=[];
                    stim_counter=0;
                    stimulus_matrix_ext=zeros(nFrames,7);
                    
                    % frame #
                    % MW/SI bitcode
                    % MW/SI bitcode
                    % trial number
                    % stimulus # (0-11)
                    % x position
                    % y position
                    % position # (1-32)
                    % position # * stim # (1-384)
                    
                    stimulus_matrix_ext(1,:)=[1 zeros(1,6)-1]; % first frame is always blank, but has no bitcode
                    glitches=0;
                    
                    for iFrame=2:nFrames
                        % step through all row of the scan image bitcodes
                        if cur_bitCode==bitCodes_scanimage(iFrame,2);
                            % do nothing
                        else
                            % if code changes copy next row of stimulus_matrix
                            cur_bitCode=bitCodes_scanimage(iFrame,2);
                            stim_counter=stim_counter+1;
                            if stim_counter>size(stimulus_matrix,1)
                                disp('ignoring end blank')
                                stim_properties=zeros(1,5)-1;
                            else
                                stim_properties=stimulus_matrix(stim_counter+glitches,3:end);
                            end
                            
                            if cur_bitCode==stim_properties(1)
                                % safe
                            else
                                glitches=glitches+1;
                                % turns out we loose frames on the SI side sometimes, I
                                % guess we should throw the corresponding MW stimulus
                                % presentations out
                                %iSess
                                %cur_bitCode
                                %stimulus_matrix(stim_counter+1,2)
                                stim_properties=stimulus_matrix(stim_counter+glitches,3:end);
                            end
                        end
                        
                        stimulus_matrix_ext(iFrame,1:7)=[iFrame cur_bitCode stim_properties];
                    end
                    %glitches
                    
                    D=diff(stimulus_matrix_ext(:,2:3),[],2);
                    
                    %%% try to cut the faulty part out
                    %stimulus_matrix_ext(diff(stimulus_matrix_ext(:,2:3),[],2)~=0,:)=[];
                    
                    if any(D)==1
                        error_row=find(D~=0,1,'first')
                        die
                        %stimulus_matrix_ext(:,2:3)
                        %error('mismatch in bitcodes, use alternative method')
                    else % skip this session
                        %% convert stim position into single number 1-32 % assumes all positions were presented
                        posMatrix=stimulus_matrix_ext(stimulus_matrix_ext(:,5)>-1,6:7);
                        positions=unique(posMatrix,'rows');
                        nPositions=size(positions,1);
                        if nPositions~=32
                            disp('Not all positions were presented, conversion will be done using static reference position matrix')
                            positions_ref=[-45.5000000000000,-19.5000000000000;-45.5000000000000,-6.50000000000000;-45.5000000000000,6.50000000000000;-45.5000000000000,19.5000000000000;-32.5000000000000,-19.5000000000000;-32.5000000000000,-6.50000000000000;-32.5000000000000,6.50000000000000;-32.5000000000000,19.5000000000000;-19.5000000000000,-19.5000000000000;-19.5000000000000,-6.50000000000000;-19.5000000000000,6.50000000000000;-19.5000000000000,19.5000000000000;-6.50000000000000,-19.5000000000000;-6.50000000000000,-6.50000000000000;-6.50000000000000,6.50000000000000;-6.50000000000000,19.5000000000000;6.50000000000000,-19.5000000000000;6.50000000000000,-6.50000000000000;6.50000000000000,6.50000000000000;6.50000000000000,19.5000000000000;19.5000000000000,-19.5000000000000;19.5000000000000,-6.50000000000000;19.5000000000000,6.50000000000000;19.5000000000000,19.5000000000000;32.5000000000000,-19.5000000000000;32.5000000000000,-6.50000000000000;32.5000000000000,6.50000000000000;32.5000000000000,19.5000000000000;45.5000000000000,-19.5000000000000;45.5000000000000,-6.50000000000000;45.5000000000000,6.50000000000000;45.5000000000000,19.5000000000000];
                            nPositions=size(positions_ref,1);
                            position_vector=zeros(nFrames,1);
                            for iPos=1:nPositions
                                pos=positions_ref(iPos,:);
                                indices=find(stimulus_matrix_ext(:,6)==pos(1)&stimulus_matrix_ext(:,7)==pos(2));
                                position_vector(indices,1)=iPos;
                            end
                        else
                            position_vector=zeros(nFrames,1);
                            for iPos=1:nPositions
                                pos=positions(iPos,:);
                                indices=find(stimulus_matrix_ext(:,6)==pos(1)&stimulus_matrix_ext(:,7)==pos(2));
                                position_vector(indices,1)=iPos;
                            end
                        end
                        
                        stimulus_matrix_ext=[stimulus_matrix_ext position_vector];
                        
                        %% convert combinations of shape and position to 1 number
                        % construct a matrix with all possible combinations of position
                        % and shape 12*32=384
                        shape_vector=stimulus_matrix_ext(stimulus_matrix_ext(:,5)>-1,5);
                        shapes=unique(shape_vector);
                        cond_matrix_full=[];
                        for iShape=1:length(shapes)
                            shape_nr=shapes(iShape);
                            cond_matrix_full=cat(1,cond_matrix_full,[repmat(shape_nr,size(positions,1),1) positions]);
                        end
                        
                        condition_vector=zeros(nFrames,1);
                        for iCond=1:size(cond_matrix_full,1)
                            condition=cond_matrix_full(iCond,:);
                            indices=find(stimulus_matrix_ext(:,5)==condition(1)&stimulus_matrix_ext(:,6)==condition(2)&stimulus_matrix_ext(:,7)==condition(3));
                            condition_vector(indices,1)=iCond;
                        end
                        %%
                        %condMatrix=stimulus_matrix_ext(stimulus_matrix_ext(:,5)>-1,[5 8]);
                        %conditions=unique(condMatrix,'rows');
                        %nConditions=length(conditions);
                        %condition_vector=zeros(nFrames,1);
                        %for iCond=1:nConditions
                        %    condition=conditions(iCond,:);
                        %    indices=find(stimulus_matrix_ext(:,5)==condition(1)&stimulus_matrix_ext(:,8)==condition(2));
                        %    condition_vector(indices,1)=iCond;
                        %end
                        stimulus_matrix_ext=[stimulus_matrix_ext condition_vector];
                        
                        file_info(matching_matrix(iSess,1)).stimulus_matrix_ext=stimulus_matrix_ext;
                        good_sessions=good_sessions+1;
                        data_sessions(good_sessions)=file_info(matching_matrix(iSess,1));
                        
                        
                    end
                catch
                    A=lasterror;
                    disp(A.message)
                    %disp('Skipped session because of mismatch in bit codes')                    
                    iSess
                end
                %size(stimulus_matrix_ext)
                
                
                if 0
                    %% Match to MWorks only based bitcodes
                    
                    %%% Create bitCode vector matching the frames by repeating each
                    %%% bitCode a given number of frames based on time between MWorks
                    %%% display update events and frame rate of SI
                    number_of_frames_per_bitCode=round((stimulus_matrix(:,2)/1e6)*frame_rate);
                    bitCode_alt=-1;
                    for iBitCode=1:length(number_of_frames_per_bitCode)
                        new=repmat(stimulus_matrix(iBitCode,3),number_of_frames_per_bitCode(iBitCode),1);
                        bitCode_alt=cat(1,bitCode_alt,new);
                    end
                    
                    %%% do check
                    if length(bitCode_alt)>nFrames
                        bitCode_alt=bitCode_alt(1:nFrames);
                        check=sum(abs(diff([stimulus_matrix_ext(:,2:3) bitCode_alt],[],2)));
                    elseif length(bitCode_alt)<nFrames
                        check=sum(abs(diff([stimulus_matrix_ext(1:length(bitCode_alt),2:3) bitCode_alt],[],2)));
                    elseif length(bitCode_alt)>size(stimulus_matrix_ext,1)
                        check=sum(abs(diff([stimulus_matrix_ext(:,2:3) bitCode_alt(1:size(stimulus_matrix_ext,1))],[],2)));
                    else
                        check=sum(abs(diff([stimulus_matrix_ext(:,2:3) bitCode_alt],[],2)));
                    end
                    check
                    if any(check)
                        error('alternative method proved to be invalid for this session...')
                    end
                end
            elseif data_type==1 % bogus flyback does not have data, so we rely soly on MWorks to identify stimulus frames
                % be careful that no other type of data falls into this category!
                %%% Match to MWorks only based bitcodes
                
                %%% Create bitCode vector matching the frames by repeating each
                %%% bitCode a given number of frames based on time between MWorks
                %%% display update events and frame rate of SI
                number_of_frames_per_bitCode=round((stimulus_matrix(:,2)/1e6)*frame_rate);
                stimulus_matrix_ext=zeros(1,7)-1;
                for iBitCode=1:length(number_of_frames_per_bitCode)
                    new=repmat([0 0 stimulus_matrix(iBitCode,3:end)],number_of_frames_per_bitCode(iBitCode),1);
                    stimulus_matrix_ext=cat(1,stimulus_matrix_ext,new);
                end
                if size(stimulus_matrix_ext,1)>nFrames
                    stimulus_matrix_ext=stimulus_matrix_ext(1:nFrames,:);
                elseif size(stimulus_matrix_ext,1)<nFrames
                    %
                end
                nFrames=size(stimulus_matrix_ext,1);
                stimulus_matrix_ext(:,1)=1:nFrames;
                
                
                %%% convert stim position into single number 1-32
                posMatrix=stimulus_matrix_ext(stimulus_matrix_ext(:,5)>-1,6:7);
                positions=unique(posMatrix,'rows');
                nPositions=size(positions,1);
                position_vector=zeros(nFrames,1);
                for iPos=1:nPositions
                    pos=positions(iPos,:);
                    indices=find(stimulus_matrix_ext(:,6)==pos(1)&stimulus_matrix_ext(:,7)==pos(2));
                    position_vector(indices,1)=iPos;
                end
                stimulus_matrix_ext=[stimulus_matrix_ext position_vector];
                
                %%% convert combinations of shape and position to 1 number
                condMatrix=stimulus_matrix_ext(stimulus_matrix_ext(:,5)>-1,[5 8]);
                conditions=unique(condMatrix,'rows');
                nConditions=length(conditions);
                condition_vector=zeros(nFrames,1);
                for iCond=1:nConditions
                    condition=conditions(iCond,:);
                    indices=find(stimulus_matrix_ext(:,5)==condition(1)&stimulus_matrix_ext(:,8)==condition(2));
                    condition_vector(indices,1)=iCond;
                end
                stimulus_matrix_ext=[stimulus_matrix_ext condition_vector];
                
                %%% Store data for later use
                file_info(matching_matrix(iSess,1)).stimulus_matrix_ext=stimulus_matrix_ext;
                good_sessions=good_sessions+1;
                data_sessions(good_sessions)=file_info(matching_matrix(iSess,1));
            end
        case 2
            disp('Parsing MWorks event data as RetinoMapping experiment')
            
            %%% Optimizations in MW protocol:
            % - include expType variable
            % - rotate bar to mark condition
            % - have frame trigger before each repeat of the sweep
            % - have more repeats
            % - have motion direction update every time so it is in
            % stimupdate event struct, so we don't have to query other
            % event codes.
                        
            %% MW: collect stimulus information
            sess_events=MW_events(MW_session_times_matched(iSess,8):MW_session_times_matched(iSess,9));
            T=double(cat(1,sess_events.time_us));
            % first stim was triggered by the end of the first SI frame so
            % appears at second frame. first frame is always free of stimulus, plus
            % it has an artifact of the delayed shutter opening => ignore
            
            nSess_events=length(sess_events);
            count=1;
            stim_counter=1;
            stimulus_matrix=[];
            last_pos_x=sess_events(1).data{2}.pos_x;
            last_pos_y=sess_events(1).data{2}.pos_y;
            for iEvent=1:nSess_events
                D=sess_events(iEvent).data; % # of session events
                timestamp=sess_events(iEvent).time_us;
                
                if length(D)==3 % stim
                    %D{2}.stim_nr=str2double(D{2}.name);
                    
                    % crude way to get information about trialtype
                    dX=D{2}.pos_x-last_pos_x;
                    dY=D{2}.pos_y-last_pos_y;
                                        
                    if D{2}.size_x<D{2}.size_y % vertical
                        if dX>0 % moving periphery to center
                            condition_nr=1;
                        else
                            condition_nr=2;
                        end
                    else % horizontal
                        if dY>0 % moving up
                            condition_nr=3;
                        else
                            condition_nr=4;
                        end
                    end
                    if dX==0&&dY==0
                        condition_nr=-1;
                    end
                    last_pos_x=D{2}.pos_x;
                    last_pos_y=D{2}.pos_y;
                    
                    %%% ideally we want orientation and direction info in
                    %%% the stim update event struct
                    %%% Rotation could be used directly to change bar and
                    %%% signal the vertical/horizontal condition. 
                    %%% Direction would then signal dX/dY
                    %%% Use isfield(D{2},'rotation') to check if variable
                    %%% is defined to distinguish between old and new format
                    
                    stimulus_matrix(count,:)=[timestamp 0 D{3}.bit_code stim_counter condition_nr D{2}.pos_x D{2}.pos_y];
                    stim_counter=stim_counter+1;
                    count=count+1;
                elseif length(D)==2 % blank, does not occur here
                    stimulus_matrix(count,:)=[timestamp 0 D{2}.bit_code -1 -1 -1 -1];
                    count=count+1;
                else
                    % empty event, can occurs at begin and end of MW session file
                end
            end
            
            if 0
                %%
                clf
                hold on
                plot(T,stimulus_matrix(:,6),'b.')
                plot(T,stimulus_matrix(:,7),'r.')
            end
            
            
            %% Now create a matrix with one row per SCIM frame and position of the stimulus per row
            % We will have to approximate this position, since the stimulus
            % moves a certain distance (frames*tempFreq=20*0.1=2degrees per scanframe)
            stimulus_matrix_ext=zeros(nFrames,7)-1;
            
            T_MW=(stimulus_matrix(:,1)-stimulus_matrix(1,1))/1e6; % construct timeline for mworks frames (60Hz)
            T_SC=(0:nFrames)/frame_rate; % construct timeline of SCIM frames
            
            for iFrame=1:nFrames
                frame_time=T_SC(iFrame);
                
                % find row in MW event matrix that is the first after the start of the SCIM frame
                row_index=find(T_MW>frame_time,1,'first');                
                
                % check time difference
                time_diff=T_MW(row_index)-frame_time;
                if time_diff>(1/frame_rate)*0.5 % blank
                    stimulus_matrix_ext(iFrame,:)=[frame_time time_diff stimulus_matrix(row_index,3:4) -1 0 0];
                else 
                    stimulus_matrix_ext(iFrame,:)=[frame_time time_diff stimulus_matrix(row_index,3:end)];
                end
                
                %stimulus_matrix_ext(iFrame,:)=[frame_time time_diff stimulus_matrix(row_index,3:end)];
            end
             
            if 0
                %%
                clf
                plot(stimulus_matrix_ext(:,1),stimulus_matrix_ext(:,2),'ro-')
            end
            
            %%
            file_info(matching_matrix(iSess,1)).stimulus_matrix_ext=stimulus_matrix_ext;
            good_sessions=good_sessions+1;
            data_sessions(good_sessions)=file_info(matching_matrix(iSess,1));
    end
end



%%% Know bad sessions
% 2015-03-04_AF11, 3 : 20150304_AF11_004
% 2015-03-05_AF03, 2 : 20150305_AF03_002
% no way we can recover this type of error with the MWorks-only method (caveat)
% 2015-04-07_AF11, 2 : 20150407_AF11_002 % some bitCodes were missed during these frames, rest of the trial, bitcodes match. Did manual correction for these

%%
if 1
    %%
        
    %%% Add the frame information only for the sessions we selected and we know are good.
    nGoodSessions=length(data_sessions);
    
    save_folder=fileparts(data_sessions(1).file_name);
    
    saveName=fullfile(save_folder,'data_analysis','session_overview.mat');
    savec(saveName)
    save(saveName,'data_sessions')
    
    
    %% Save individual sessions: will destroy existing files e.g. motion correction, ROI_definitions
    t0=clock;
    for iSess=1:nGoodSessions
        session_data=data_sessions(iSess);
        [save_folder, save_name]=fileparts(session_data.file_name);
        saveName=fullfile(save_folder,'data_analysis',[save_name '.mat']);
        save(saveName,'session_data')
        progress(iSess,nGoodSessions,t0)
    end
    disp('All Done!!!')
else
    disp('Not saving for now, debugging')
end


