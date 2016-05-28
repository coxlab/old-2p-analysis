clear all
clc


%png_name='/Volumes/TeraLacie/eyeTracker/Calibration_2015-08-19/left_137.png';
%png_name='/Volumes/TeraLacie/eyeTracker/Calibration_2015-08-19/right_137.png';

%png_name='/Volumes/TeraLacie/eyeTracker/2015-08-19_AH03/.png';

%im=double(imread(png_name));

%imshow(im,[])

load_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-19_AH03/data_analysis/2015-08-19_AH03_001.mat';
load(load_name,'session_data')
session_data

im_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/eyeTracker/2015-08-19_AH03';
tic
files=scandir(im_folder,'png');
nFiles=length(files);
toc

%%
tic

run(nFiles)=struct('run_ID','','run_nr',[],'frame_nr_est',[],'s1',[],'s2',[],'ts',[],'frame_nr_real',[],'bitCode',[]);
cur_run_id='';
cur_run_nr=0;
toc


tic
for iFile=1:nFiles
     %im_name=fullfile(im_folder,files(iFile).name);
     % parse filename
     [f,file_name]=fileparts(files(iFile).name);
     file_parts=strsplit(file_name,'_');
     
     run_ID=file_parts{2};
     if strcmpi(run_ID,cur_run_id)==0
         cur_run_id=run_ID;
         cur_run_nr=cur_run_nr+1;
     end
     
     run(iFile).run_ID=run_ID;
     run(iFile).run_nr=cur_run_nr;
     run(iFile).frame_nr_est=str2double(file_parts{3});
     run(iFile).s1=strcmpi(file_parts{5},'1');
     run(iFile).s2=strcmpi(file_parts{7},'1');
     run(iFile).ts=str2double(file_parts{9});
     run(iFile).frame_nr_real=str2double(file_parts{4});
     run(iFile).bitCode=sum([run(iFile).s1 run(iFile).s2].*[1 2]);
end
toc
%%

run_nr_vector=cat(1,run.run_nr);
frame_nr_vector=cat(1,run.frame_nr_est);
tabulate(run_nr_vector)

%%
sel=run_nr_vector==1;
V=frame_nr_vector(sel);
plot(diff(sort(V)))

%%
V_sorted=sort(V);
diff(V_sorted(1:100))

%% build matrix with all relevant info
run_nr_vector=cat(1,run.run_nr);
frame_nr_vector=cat(1,run.frame_nr_est);
s1_vector=cat(1,run.s1);
s2_vector=cat(1,run.s2);
bitCode_vector=cat(1,run.bitCode);
%%
M=[run_nr_vector frame_nr_vector s1_vector s2_vector bitCode_vector];
M=sortrows(M,[1 2]);
T=(1:size(M,1))'/block_size_est;
M=[T M];

for run_nr=1%:6
    D=M(M(:,2)==run_nr,:);
    
    subplot(3,2,run_nr)
    plot(D(:,1),D(:,end))
end 

ET_bitCodes_raw=D(:,end);

%%
block_size_est=15;
nBlocks=floor(length(ET_bitCodes_raw)/block_size_est)-1;
offset=4;
ET_bitCodes=mode(reshape(ET_bitCodes_raw(offset+1:offset+nBlocks*block_size_est),block_size_est,[]))';

%% find switches and times
%plot(ET_bitCodes)


%%
load_name='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-19_AH03/data_analysis/2015-08-19_AH03_001.mat';
load(load_name,'session_data')
MW_bitCodes=session_data.bitCodes.MWorks_bitCodes(:,2);

%plot(mod(MW_bitCodes,4))

%% expand to two bitCodes in case of blank trial 2s interval...
MW_bitCodes_2bit=mod(MW_bitCodes,4);


T=session_data.bitCodes.MWorks_bitCodes(:,1);
%T=T-T(1);
T=T-T(4);
durations=[0 ; diff(session_data.bitCodes.MWorks_bitCodes(:,1))];

M_MW=[T MW_bitCodes_2bit durations];

%M_MW([1:3 end],:)=[];

subplot(211)
bar(D(:,1),D(:,end)+1)
axis([450 500 0 4])
subplot(212)
bar(M_MW(:,1),M_MW(:,2)+1)
axis([7740 7800 0 4])


%% Parse up mworks bitcodes to have equally spaced samples
first_stim=find(between(durations,[.95 1.05]*1),1,'first');
last_blank=find(between(durations,[.95 1.05]*2),1,'last');

M_sel=M_MW(first_stim:last_blank,:);

M_reshape=reshape(M_sel(:,2),2,[]);

M_reshape=cat(1,M_reshape,M_reshape(2,:));

M_1hz=M_reshape(:);

%%

%CC=normxcorr2(ET_bitCodes(:),M_1hz(:));

[cc_value,offset]=im_align(ET_bitCodes,M_1hz);


matched=[ET_bitCodes M_1hz(abs(offset(2))+1:abs(offset(2))+length(ET_bitCodes))];
matched(1:20,:)

%unique({run.run_ID},'stable')'
%
% Value    Count   Percent
%       1    22175     20.64%
%       2    17097     15.91%
%       3    16951     15.77%
%       4    16504     15.36%
%       5    17505     16.29%
%       6    17231     16.03%

% In order newest to oldest
%     '548164372a6969bc1bc26937ff1facf7'
%     '571dc0ce854cc1f9dd5025a64feb132a'
%     '60a0b71f20090ee84cbf209934a2e809'
%     '84e61962f39cd2e94297b308ce611b2c'
%     'a99bb417b58ed3c567c5241bc9d2cf7c'
%     'f1733442681fc30cf07fb92e5956b793'

