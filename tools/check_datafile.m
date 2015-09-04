clear all
clc

switch 2
    case 1
        datafile_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-16_AF11';
        file_name=fullfile(datafile_folder,'20150416_AF11_001.tif')
    case 2
        datafile_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-07-17_AG02_awake';
        file_name=fullfile(datafile_folder,'2015-07-17_AG02_003.tif')
end

info=imfinfo(file_name);

%%% Get basic info about movie
nFrames=length(info);
W=info(1).Width;
H=info(1).Height;
scim_info=info(1).ImageDescription;
scinfo=strsplit(scim_info,char(13));
scim_info=strjoin(scinfo,[';' char(13)]);
eval([scim_info  ';']); % revive scan image variables
frame_rate=state.acq.frameRate;


N=[];
%frame_info=struct;
for iFrame=1:nFrames %min([500 nFrames])
    % read out last line of the frame, here we store frame specific
    % info: bitCodes, position, laser power
    flyback_line=double(imread(file_name,iFrame,'info',info,'PixelRegion',{[H H],[1 W]}));
    
    [a,N]=parse_flyback_line(flyback_line,N);
    frame_info(iFrame)=a;
    
end

%% sanity checks
A=cat(1,frame_info.nBitCodes);
if all(ismember(A,mode(A)))
    disp('Passed!')
else
    error('variable number of bitCodes detected')
end


%%
A=cat(1,frame_info.timestamp);
%A=A-A(1);
D=diff(A);

%ttest_report(D,1/frame_rate)
[mean(D) std(D) 1/frame_rate]

Z=zscore(D);
plot(Z)
if sum(any(Z>4))>2
    error('Timing is off')
else
    disp('Passed interval test')
end
if 0
    %%
    hist(D,length(A)/20)
    hold on
    plot(1./[frame_rate frame_rate],[0 100 ],'r')
    hold off
end
%%
%%% Unwrap full bitCode vector
A=cat(1,frame_info.bitCode_vector);
A=medfilt1(A);

breaks=find([-1 ; diff(A)]);
sel=diff([breaks;0])<149;
if length(sel)==length(breaks)
    breaks(sel)=[];
else
    error('messed up...')
end

bitCode_time=(1/frame_rate)/N;
T=((0:length(A)-1))'*bitCode_time;

scanimage_bitCodes=[T(breaks) A(breaks)];
D=diff(scanimage_bitCodes(:,1));

hist(D,length(D)/10)
size(scanimage_bitCodes)

if 0
    %%
    plot(T,A)
    hold on
    plot(T(breaks),ones(size(T(breaks))),'m*')
    hold off
end

%%
tic
mwk_file_name=fullfile(datafile_folder,'2015-07-17_AG02.mwk');
disp('Reading MWK file...')
A=getCodecs(mwk_file_name);
event_codec=A.codec;

%%% Get stim update events
tag_name='#stimDisplayUpdate';
[MW_events,nEvents,event_code]=get_events_by_name(mwk_file_name,tag_name,event_codec);

MW_bitCodes=zeros(nEvents,3);
for iEvent=1:nEvents
    event=MW_events(iEvent);
    register=1;
    if length(event.data)==3
        bitCode=event.data{3}.bit_code;
    elseif length(event.data)==2
        bitCode=event.data{2}.bit_code;
    else
        % ignore
        register=0;
    end
    if register==1
        MW_bitCodes(iEvent,:)=[iEvent double(event.time_us)/1e6 double(bitCode)];
    end
end

MW_bitCodes(:,[1 3])

fprintf('Loading MWorks events took %3.2f seconds.\n',toc)

%% Matching
A=scanimage_bitCodes(:,2);
B=MW_bitCodes(:,3);

CC=normxcorr2(A,B);
[maxVal,loc]=max(CC);
offset=loc-length(A)+1;
check=[A B(offset:offset+length(A)-1) diff([A B(offset:offset+length(A)-1)],1,2)];
if any(check(:,3)~=0)
    disp('Solution is not pure...')
    check(check(:,3)~=0,:)
end

plot(CC)
title(sprintf('max corr=%3.2f',maxVal))
xlabel(offset)
axis([length(A) length(B) -1 1])

if 0
    %% Validate offset
    diff([MW_events([offset-2 offset-1 offset offset+1 offset+2 offset+3]).time_us])
end



