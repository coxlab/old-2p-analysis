clear all
clc

% load ROIs from different files and check their alignment.
% additional info to be used
% - FOV position
% - coregistration of MIP in both, can be different size

header_script

folder_selector=1;
cluster_nr=1;

%%% Select datafiles
use_gui=0;
switch folder_selector
    case 0
        %data_folder=uigetdir(data_root);
        cd(data_root)
        switch 2
            case 1 % handpick session recorded at same FOV and using same stimulation protocol
                [filenames, data_folder]=uigetfile('.mat','Pick data file','MultiSelect', 'on');
                use_gui=1;
            case 2 % select folder and use FOV_matching variable in session_overview to select session that belong together
                % this is cool, but we need another matrix that specifies
                % experiment type so we can focus on same FOV and same exp
                %data_folder=uigetdir(data_root);
                
                loadName=fullfile(data_folder,'data_analysis','session_overview.mat');
                load(loadName,'data_sessions','FOV_matching')
                
                nClusters=max(FOV_matching.clusters(:));
                if cluster_nr>nClusters
                    cluster_nr=nClusters;
                    disp('Reached max number of clusters');
                end
                [r,c]=find(FOV_matching.clusters==cluster_nr);
                selected_sessions=unique(r);
                
                filenames=cell(0,1);
                for iSess=1:length(selected_sessions)
                    [f,fn,ext]=fileparts(data_sessions(selected_sessions(iSess)).file_name);
                    filenames{iSess}=[fn '.mat'];
                end
                data_folder=fullfile(data_folder,'data_analysis');
                use_gui=1;
        end
        %loadName=fullfile(data_folder,filename);
        %session_nr_vector=[1];
        %loadName_format='20150407_jrat3_%03d.mat';
        
    case 1
        data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-05_AF11/data_analysis/';
        
        loadName_format='20150305_AF11_%03d.mat';
        cluster_nr=3;
        switch cluster_nr
            case 1
                session_nr_vector=[1 15 17];
            case 2
                session_nr_vector=[2 3 4]; % 18
            case 3
                session_nr_vector=[8 9]; % 12Hz
            case 4 % cell 8 shows a sign of RF
                session_nr_vector=[10 11 12]; % 12Hz
        end
    case 2
        data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-03-04_AF11/data_analysis/';
        loadName_format='20150304_AF11_%03d.mat';
        cluster_nr=1;
        session_nr_vector=[2 3 5];
    case 3
        data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-07_AF11/data_analysis/';
        loadName_format='20150407_AF11_%03d.mat';
        cluster_nr=1;
        switch cluster_nr
            case 1
                session_nr_vector=[1 2];
            case 2
                session_nr_vector=[3 4 5];
            case 3
                session_nr_vector=[6 7];
        end
end



%%
%ROI_all=struct;
if use_gui==1
    nSessions=length(filenames);
    for iSess=1:nSessions
        loadName=fullfile(data_folder,filenames{iSess});
        if exist(loadName,'file')
            load(loadName,'session_data')
            if isfield(session_data,'options')
                % remove old options field
                session_data=rmfield(session_data,'options');
            end
            ROI_all(iSess)=session_data;
        end
    end
else
    nSessions=length(session_nr_vector);
    for iSess=1:nSessions
        session_nr=session_nr_vector(iSess);
        loadName=fullfile(data_folder,sprintf(loadName_format,session_nr));
        
        if exist(loadName,'file')
            load(loadName,'session_data')
            ROI_all(iSess)=session_data;
        end
    end
end


%%
%%% given we have copied ROIs, we can be sure that same numbers indicate
%%% same cells. Deleting ROI and creating new ones results in new unique
%%% numbers.

% get all unique numbers, find in all repetitions of the experiment
unique_cell_numbers=[];
for iSess=1:nSessions
    cell_numbers=cat(1,ROI_all(iSess).ROI_definitions.ROI_nr);
    %all_cell_numbers=cat(1,all_cell_numbers,cell_numbers);
    if isempty(unique_cell_numbers)
        unique_cell_numbers=cell_numbers;
    else
        unique_cell_numbers=intersect(unique_cell_numbers,cell_numbers);
    end
end

nROI=length(unique_cell_numbers);

%%% Glue session together
STIM_ALL=[];
RESP_ALL=[];
RESP_ALL_oopsi=[];
data_matrix_all=[];
data_matrix_base_all=[];
data_matrix_all_oopsi=[];
data_matrix_base_all_oopsi=[];

for iSess=1:nSessions
    cell_numbers=cat(1,ROI_all(iSess).ROI_definitions.ROI_nr);
    cell_selection=ismember(cell_numbers,unique_cell_numbers);
    time_selection=2:size(ROI_all(iSess).stimulus_matrix_ext,1);
    
    STIM=ROI_all(iSess).stimulus_matrix_ext(time_selection,:);
    RESP=ROI_all(iSess).activity_matrix(time_selection,cell_selection);
    %RESP=ROI_all(iSess).spike_matrix(time_selection,cell_selection);
    RESP_oopsi=ROI_all(iSess).spike_matrix(time_selection,cell_selection);
    %[size(STIM) size(RESP)]
    STIM_ALL=cat(1,STIM_ALL,STIM);
    RESP_ALL=cat(1,RESP_ALL,RESP);
    RESP_ALL_oopsi=cat(1,RESP_ALL_oopsi,RESP_oopsi);
    
    % chop up RESP into trial responses using STIM
    trial_matrix=STIM(STIM(:,4)>0,[1 4]);
    
    % find frames where new trial begins
    new_trial_frames=[2; trial_matrix([0;diff(trial_matrix(:,2))]==1,1)]; % new trial at start of exp
    
    nTrials=max(trial_matrix(:,2));
    
    data_matrix=zeros(nTrials-2,10+nROI);
    data_matrix_base=data_matrix;
    data_matrix_oopsi=data_matrix;
    data_matrix_base_oopsi=data_matrix;
    for iTrial=2:nTrials-1  % ignore first trial, might be infected with laser not being on, ignore last trial, because it could cut off
        switch 4
            case 1 % only stim frames
                sel=trial_matrix(:,2)==iTrial;
            case 2 % stim+blank after
                sel=new_trial_frames(iTrial):new_trial_frames(iTrial+1)-1;
            case 3
                stim_window=round((new_trial_frames(iTrial+1)-new_trial_frames(iTrial))*0.75);
                sel=new_trial_frames(iTrial):new_trial_frames(iTrial)+stim_window;
            case 4
                %stim_window=round((new_trial_frames(iTrial+1)-new_trial_frames(iTrial))*0.75);
                stim_window=round(round((new_trial_frames(iTrial+1)-new_trial_frames(iTrial))*[0.25 0.75]));
                base_window=round(round((new_trial_frames(iTrial+1)-new_trial_frames(iTrial))*[-0.25 0.25]));
                sel=new_trial_frames(iTrial)+stim_window(1):new_trial_frames(iTrial)+stim_window(2);
                sel_base=new_trial_frames(iTrial)+base_window(1):new_trial_frames(iTrial)+base_window(2);
        end
        
        stim_properties=STIM(STIM(:,1)==new_trial_frames(iTrial),[1 3:end]);
        
        %%% Regular signal
        trial_data=RESP(sel,:); % all frames collected during a given trial
        trial_data_base=RESP(sel_base,:); % all frames collected during a given trial
        data_matrix(iTrial-1,:)=[iSess iTrial stim_properties mean(trial_data,1)];
        data_matrix_base(iTrial-1,:)=[iSess iTrial stim_properties mean(trial_data_base,1)];
        
        %%% Non negative deconvoluted signal
        trial_data=RESP_oopsi(sel,:); % all frames collected during a given trial
        trial_data_base=RESP_oopsi(sel_base,:); % all frames collected during a given trial
        data_matrix_oopsi(iTrial-1,:)=[iSess iTrial stim_properties mean(trial_data,1)];
        data_matrix_base_oopsi(iTrial-1,:)=[iSess iTrial stim_properties mean(trial_data_base,1)];
    end
    
    data_matrix_all=cat(1,data_matrix_all,data_matrix);
    data_matrix_base_all=cat(1,data_matrix_base_all,data_matrix_base);
    data_matrix_all_oopsi=cat(1,data_matrix_all_oopsi,data_matrix_oopsi);
    data_matrix_base_all_oopsi=cat(1,data_matrix_base_all_oopsi,data_matrix_base_oopsi);
end

%size(STIM_ALL)
%size(RESP_ALL)



%%
% 20150505: do spike detection on traces and sort ROIs by descending number of spikes. Did we ever implement this function?
TH=5;
for iROI=1:nROI
    sel=RESP_ALL(:,iROI)>TH;
    start_points=find(diff(sel)==1);
    nSpikes_vector(iROI)=length(start_points);
end
[sorted,order]=sort(nSpikes_vector,'descend')
cell_numbers_ranked=unique_cell_numbers(order);













%%
%folder_selector=1
% cluster 1: 12? 26shape 35shape? 43?
% cluster 2: 1? 11pos 16shape
% cluster 3: 1? 4? 6 7 8 11axon 13!axon 14axon 16?axon 17?axon 19?axon 22? 23axon 24!shape&pos
% cluster 4 s: 3pos&shape 5? 6? 7? 8 9?axon&shape 10 11 13shape 14pos&shape 16?shape 17 20?axon 22?axon 23? 28shape 33pos 34?pos&shape

%folder_selector=2
% cluster 1: c2shape? c16shape? c19?pos c21pos c23pos? 32?pos
plot_it=1;
ROI_sel=9
use_NND=0;

ROI_nr=cell_numbers_ranked(ROI_sel);
if ismember(ROI_nr,unique_cell_numbers)
    iROI=find(unique_cell_numbers==ROI_nr); % Cluster3 11 13 24
else
    error('no such ROI nr')
end

%ROI_nr=session_data.ROI_definitions(iROI).ROI_nr

%%% look at response per stimulus condition
% average all fluo responses for stim and blank after (3s)

nPositions=32;
resp_data=zeros(nPositions,nROI,3);
base_data=resp_data;
Y=[];groups=[];
for iPos=1:nPositions
    sel=data_matrix_all(:,9)==iPos;
    switch use_NND
        case 0
            resp_matrix=data_matrix_all(sel,11:end);
            base_matrix=data_matrix_base_all(sel,11:end);
        case 1
            resp_matrix=data_matrix_all_oopsi(sel,11:end);
            base_matrix=data_matrix_base_all_oopsi(sel,11:end);
    end
    
    %nTrials=size(resp_matrix,1);
    resp_data(iPos,:,:)=cat(3,mean(resp_matrix,1),std(resp_matrix,[],1),ste(resp_matrix,1));
    base_data(iPos,:,:)=cat(3,mean(base_matrix,1),std(base_matrix,[],1),ste(base_matrix,1));
    net_data=cat(3,mean(resp_matrix,1)-mean(base_matrix,1),std(resp_matrix,[],1),ste(resp_matrix,1));
    
    % build anova matrix
    Y=cat(1,Y,resp_matrix);
    groups=cat(1,groups,repmat(iPos,size(resp_matrix,1),1));
end
[p_pos,t,stats]=anovan(Y(:,iROI),groups,'display','off');

if plot_it==1
    figure(2)
    subplot(221)
    bar(resp_data(:,iROI,1))
    hold on
    %errorbar(resp_data(:,iROI,1),resp_data(:,iROI,2),'r.')
    errorbar(resp_data(:,iROI,1),resp_data(:,iROI,2),'r.')
    hold off
    set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
    title('Screen Position')
    xlabel(sprintf('Anova p=%3.4f',p_pos))
    axis([1 nPositions -0.5 max(resp_data(:,iROI,1))*2])
    
    subplot(222)
    RF=flipud(reshape(resp_data(:,iROI,1),4,8));
    max_val=max(abs([min(RF(:)) max(RF(:))]));
    imagesc(RF)
    axis equal
    axis tight
    
    %set(gca,'CLim',[-max_val max_val],'xTickLabel',[],'yTickLabel',[])
    set(gca,'CLim',[-.5 .5],'xTickLabel',[],'yTickLabel',[])
    title('''RF''')
    set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
    colormap(parula)
    colorbar
end

% Examine responses
stim_vector=STIM_ALL(:,4)>0;
resp_vector=RESP_ALL(:,iROI);
resp_vector_oopsi=RESP_ALL_oopsi(:,iROI);

nFrames=size(stim_vector,1);

%%% Construct time scale
frame_rate=session_data.data(5);
T=((1:nFrames)-1)/frame_rate;
total_duration=max(T);


[sorted,order]=sort(resp_data(:,iROI,1),'descend');
nBest=2;
best_conditions=ismember(STIM_ALL(:,8),order(1:nBest));
best_condition=ismember(STIM_ALL(:,8),order(1));
sorted_positions=order;
if plot_it==1
    
    
    
    subplot(2,2,[3 4])
    max_val=max([30 max(resp_vector)*1.2]);
    bar(T,spherify(stim_vector,2)*max_val,'barWidth',1,'FaceColor',[1 1 1]*.9,'EdgeColor',[1 1 1]*.9)
    hold on
    bar(T,spherify(best_conditions,2)*-1,'barWidth',1,'FaceColor',[1 0 0]*.5,'EdgeColor',[1 0 0]*.5)
    bar(T,spherify(best_condition,2)*-1,'barWidth',1,'FaceColor',[1 0 0]*.9,'EdgeColor',[1 0 0]*.9)
    plot(T,resp_vector,'color','k')
    plot(T,resp_vector_oopsi,'color','g')
    hold off
    axis([0 total_duration -max_val*.25 max_val])
    set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
    title(sprintf('Cell #%d - #spikes: %d',[ROI_nr nSpikes_vector(iROI)]))
    xlabel('Time (seconds)')
    ylabel('\DeltaF/F')
end

%%% double check
% sel=data_matrix_all(:,8)==order(2);
% resp_matrix=data_matrix_all(sel,10:end);
% resp_matrix(:,iROI)
% nTrials=size(resp_matrix,1);
% new=cat(3,mean(resp_matrix,1)-1,std(resp_matrix,[],1),ste(resp_matrix,1));


%%
timepoint_vector=1:9;
nTimepoints=length(timepoint_vector);
condition_vector=STIM_ALL(:,8);
switch 1
    case 1
        %% find kernel estimate using different timepoints of all conditions
        % find all condition start frames
        start_vector=find(diff(condition_vector)>0); % general case
        sel=(start_vector+max(timepoint_vector))>nFrames|(start_vector+min(timepoint_vector))<0;
        start_vector(sel)=[];
        
        design_matrix=ones(size(resp_vector));
        for iTimepoint=1:nTimepoints
            P=zeros(size(resp_vector));
            P(start_vector+timepoint_vector(iTimepoint))=1;
            design_matrix=cat(2,design_matrix,P);
        end
        beta=mvregress(design_matrix,resp_vector);
        kernel_est=beta(2:end);
        %kernel_est(kernel_est<0)=0;
        
    case 2
        %% find kernel estimate using different timepoints of all conditions
        % find all condition start frames
        start_vector=find(diff(condition_vector)>0); % general case
        sel=(start_vector+max(timepoint_vector))>nFrames|(start_vector+min(timepoint_vector))<0;
        start_vector(sel)=[];
        
        design_matrix=ones(size(resp_vector));
        for iTimepoint=1:nTimepoints
            P=zeros(size(resp_vector));
            P(start_vector+timepoint_vector(iTimepoint))=1;
            design_matrix=cat(2,design_matrix,P);
        end
        beta=mvregress(design_matrix,resp_vector);
        kernel_est=beta(2:end);
        %kernel_est(kernel_est<0)=0;
        
    case 3
        %% find kernel estimate using different timepoints, for each condition
        %condition_vector=STIM_ALL(:,8);
        kernel_est_all=zeros(nTimepoints,nPositions);
        for iStim=1:nPositions
            P=condition_vector==iStim;
            
            %%% estimate kernel for each condition
            start_vector=find(diff(condition_vector)==iStim);
            sel=(start_vector+max(timepoint_vector))>nFrames|(start_vector+min(timepoint_vector))<0;
            start_vector(sel)=[];
            
            design_matrix_stim=ones(size(resp_vector));
            for iTimepoint=1:nTimepoints
                P=zeros(size(resp_vector));
                P(start_vector+timepoint_vector(iTimepoint))=1;
                design_matrix_stim=cat(2,design_matrix_stim,P);
            end
            beta=mvregress(design_matrix_stim,resp_vector);
            kernel_est_all(:,iStim)=beta(2:end);
        end
        
        [m,w]=max(max(kernel_est_all));
        kernel_est=kernel_est_all(:,w);
end
figure(4)
bar(kernel_est)


%%

%%% try out multiple regression
%kernel_est=resp_vector(365:390);
%kernel_est=[34.7974745645518;29.1844181233045;24.5764182663738;22.8407164563818;18.1783359501114;16.4113462293559;13.4066753860348;14.1086533821770;10.1052434699760;10.0838630155124;8.25588908292035;6.43065557040027;6.04692099511076;5.16536561771791;4.76037094647263;6.56684814422224;3.81085714655538;2.92912801412736;3.03125673618300;1.93144429292142;2.00843847510978;2.23100586840692;0.639220487381043;0.459113580149153;1.01104602349656;0.234196624416639];
design_matrix=ones(size(resp_vector));
%condition_vector=STIM_ALL(:,8);
for iStim=1:32
    P=condition_vector==iStim;
    
    %%% use average kernel
    P_conv=convn(P,kernel_est);
    P_conv(end-nTimepoints+2:end)=[]; % trim convoluted vector
    
    %%% Do multiple regression
    design_matrix=cat(2,design_matrix,P_conv);
end
% plot(P)
% hold on
% plot(P_conv,'r')
% hold off

beta=mvregress(design_matrix,resp_vector);

RF_regr=flipud(reshape(beta(2:end),4,8));


CC=corr(resp_vector,design_matrix);
RF_corr=flipud(reshape(CC(2:end),4,8));
%RF_max_kernel=flipud(reshape(max(kernel_est_all),4,8));

[m,w]=max(CC(:));

%%
figure(3)
subplot(2,3,[1 3])
plot(T,resp_vector)
hold on
plot(T,design_matrix(:,w)/max(design_matrix(:,w))*max(resp_vector),'r')
plot(T,resp_vector_oopsi,'g')
hold off
title(corr(resp_vector,design_matrix(:,w)))
max_val=max([1 max(resp_vector)*1.2]);
%axis([0 total_duration -1 1])
axis([0 total_duration -max_val/4 max_val])
set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})

subplot(234)
max_val=max(abs([min(RF_regr(:)) max(RF_regr(:))]));
imagesc(RF_regr)
axis equal
axis tight

set(gca,'CLim',[-max_val max_val],'xTickLabel',[],'yTickLabel',[])
title('''RF multiple regression''')
colormap parula
colorbar

subplot(235)
max_val=max(abs([min(RF_corr(:)) max(RF_corr(:))]));
imagesc(RF_corr)
axis equal
axis tight

set(gca,'CLim',[-max_val max_val],'xTickLabel',[],'yTickLabel',[])
title('RF correlation')
colormap parula
colorbar

subplot(236)
cla
% max_val=max(abs([min(RF_max_kernel(:)) max(RF_max_kernel(:))]));
% imagesc(RF_max_kernel)
% axis equal
% axis tight
% set(gca,'CLim',[-max_val max_val],'xTickLabel',[],'yTickLabel',[])
% title('RF max kernel')
% colormap parula
% colorbar


%%


%% Effect of Shape Stimulus

nShapes=12;
resp_data=zeros(nShapes,nROI,3);
Y=[];groups=[];
for iShape=1:nShapes
    sel=data_matrix_all(:,6)==iShape-1;
    resp_matrix=data_matrix_all(sel,11:end);
    nTrials=size(resp_matrix,1);
    resp_data(iShape,:,:)=cat(3,mean(resp_matrix,1),std(resp_matrix,[],1),ste(resp_matrix,1));
    
    % build anova matrix
    Y=cat(1,Y,resp_matrix);
    groups=cat(1,groups,repmat(iShape,size(resp_matrix,1),1));
end
size(resp_data)
[p_shape,t,stats]=anovan(Y(:,iROI),groups,'display','off');


%%%
figure(3)

subplot(2,2,[1 2])
bar(resp_data(:,iROI,1))
hold on
errorbar(resp_data(:,iROI,1),resp_data(:,iROI,2),'r.')
errorbar(resp_data(:,iROI,1),resp_data(:,iROI,3),'g.')
hold off
set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
title('Shape')
xlabel(sprintf('Anova p=%3.4f',p_shape))

% Examine responses
stim_vector=STIM_ALL(:,4)>0;
resp_vector=RESP_ALL(:,iROI);
nFrames=size(stim_vector,1);

[sorted,order]=sort(resp_data(:,iROI,1),'descend');
nBest=1;
best_conditions=ismember(STIM_ALL(:,5),order(1:nBest)-1);
best_condition=ismember(STIM_ALL(:,5),order(1)-1);
order(1:nBest)

subplot(2,2,[3 4])
%max_val=max([4 max(resp_vector)]);
max_val=max([30 max(resp_vector)*1.2]);
bar(T,spherify(stim_vector,2)*max_val,'barWidth',1,'FaceColor',[1 1 1]*.9,'EdgeColor',[1 1 1]*.9)
hold on
bar(T,spherify(best_conditions,2)*-1,'barWidth',1,'FaceColor',[1 0 0]*.5,'EdgeColor',[1 0 0]*.5)
bar(T,spherify(best_condition,2)*-1,'barWidth',1,'FaceColor',[1 0 0]*.9,'EdgeColor',[1 0 0]*.9)

plot(T,resp_vector,'color',[0 0 0]*.7)
axis([1 total_duration -max_val*.25 max_val])
set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
hold off
title(sprintf('Cell #%d',ROI_nr))
%xlabel('Frame number')
xlabel('Time (seconds)')
ylabel('\DeltaF/F')


%%% Response to shape for different position
figure(4)
if 0
    %%
    condition_vector=data_matrix_all(:,9);
    conditions=unique(condition_vector);
    nConditions=384;
    
    resp_data=zeros(nConditions,nROI,2);
    Y=[];groups=[];
    for iCond=1:nConditions
        sel=data_matrix_all(:,9)==iCond;
        if sum(sel)>0
            resp_matrix=data_matrix_all(sel,10:end);
            resp_data(iCond,:,:)=cat(3,mean(resp_matrix,1),std(resp_matrix,[],1));
        end
    end
end


%if p_pos>.05
nBest=1;
% next best thing: check stimulus response on the best positions,
% cutting out the trials where stim is outside 'RF'
resp_data=zeros(nShapes,nROI,2);
for iShape=1:nShapes
    sel=data_matrix_all(:,6)==iShape-1&ismember(data_matrix_all(:,9),sorted_positions(1:nBest));
    if sum(sel)>0
        resp_matrix=data_matrix_all(sel,11:end);
        resp_data(iShape,:,:)=cat(3,mean(resp_matrix,1),std(resp_matrix,[],1));
    end
end
bar(resp_data(:,iROI,1))
hold on
errorbar(resp_data(:,iROI,1),resp_data(:,iROI,2),'r.')
hold off
%else
%    clf
%end


%% now, pick to shapes and plot average dF/F for both, using 12Hz data
if cluster_nr==3
    iROI=13
    figure(5)
    switch iROI
        case 11
            shape_pair=[9 11]; % for 3-11
            nPositions=3;
            max_val=4;
        case 13
            shape_pair=[7 9 10]; % for 3-13
            nPositions=2;
            max_val=20;
            nROI
        case 24
            shape_pair=[6 10]; % for 3-13
            nPositions=3;
            max_val=4;
        otherwise
            die
    end
    nFrames=81;
    pre_frames=20;
    selected_positions=sorted_positions(1:nPositions);
    selected_positions
    
    
    % preferred shape
    sel=data_matrix_all(:,6)==shape_pair(3)-1&ismember(data_matrix_all(:,9),selected_positions);
    trial_nrs=data_matrix_all(sel,1:3);
    P=[];
    for iFrame=1:size(trial_nrs,1)
        start_frame=find(STIM_ALL(:,1)==trial_nrs(iFrame,3)&STIM_ALL(:,4)==trial_nrs(iFrame,2));
        P(iFrame,:)=RESP_ALL(start_frame-pre_frames:start_frame+nFrames-pre_frames-1,iROI);
    end
    
    % neutral shape
    sel=data_matrix_all(:,6)==shape_pair(2)-1&ismember(data_matrix_all(:,9),selected_positions);
    trial_nrs=data_matrix_all(sel,1:3);
    N=[];
    for iFrame=1:size(trial_nrs,1)
        start_frame=find(STIM_ALL(:,1)==trial_nrs(iFrame,3)&STIM_ALL(:,4)==trial_nrs(iFrame,2));
        N(iFrame,:)=RESP_ALL(start_frame-pre_frames:start_frame+nFrames-pre_frames-1,iROI);
    end
    
    
    % non-preferred shape
    sel=data_matrix_all(:,6)==shape_pair(1)-1&ismember(data_matrix_all(:,9),selected_positions);
    trial_nrs=data_matrix_all(sel,1:3);
    NP=[];
    for iFrame=1:size(trial_nrs,1)
        start_frame=find(STIM_ALL(:,1)==trial_nrs(iFrame,3)&STIM_ALL(:,4)==trial_nrs(iFrame,2));
        NP(iFrame,:)=RESP_ALL(start_frame-pre_frames:start_frame+nFrames-pre_frames-1,iROI);
    end
    X=(-pre_frames:nFrames-pre_frames-1)/frame_rate;
    %p(1)=errorbar(X,mean(P,1),ste(P,1),'k.-');
    %
    M=[];
    M(:,1)=mean(P,1);
    M(:,2)=mean(N,1);
    M(:,3)=mean(NP,1);
    switch 2
        case 1
            for i=1:3
                M(:,i)=medfilt1(M(:,i),6);
            end
        case 2
            G=bellCurve(1,0,frame_rate/4,round(frame_rate))';
            G=G./sum(G(:));
            M=convn(M,G);
            M=M(size(G,1)/2:end-size(G,1)/2,:);
        case 3
            %%
            F=M(:,1);
            G=bellCurve(1,0,length(F)/2,length(F))';
            %G=G-mean(G);
            DC=min(F(:));
            F=F-DC;
            FFT=fft(F);
            %DC=FFT(1);
            FFT=FFT-DC;
            F_fft=-real(ifft(fftshift(fftshift(FFT).*G)));
            M(:,1)=F_fft;
    end
    
    %%% plot shit
    p(1,:)=shadedErrorBar(X,M(:,1),ste(P,1),'r-');
    hold on
    p(2,:)=shadedErrorBar(X,M(:,2),ste(N,1),'g-');
    p(3,:)=shadedErrorBar(X,M(:,3),ste(NP,1),'b-');
    line([0 0],[0 max_val],'color','r')
    line([2 2],[0 max_val],'color','k')
    line([3 3],[0 max_val],'color','r')
    hold off
    axis([-1 5 0 max_val])
    box off
    title('Average trace for good and bad shape')
    xlabel('Time (seconds)')
    ylabel('\DeltaF/F')
    %ylabel('Average Fluorescence')
    legend(cat(1,p.mainLine),{sprintf('Preferred (N=%d)',size(P,1)),sprintf('Neutral (N=%d)',size(N,1)),sprintf('Non-Preferred (N=%d)',size(NP,1))})
    
    print(gcf,'traces.eps','-depsc')
end















%%
if 0
    %% check whether field overlap
    for iSess=1:length(session_nr_vector)
        if isfield(ROI_all(iSess).FOV_info,'center')
            ROI_all(iSess).FOV_info
        end
    end
    
    
    %% calculate all pairwise shifts
    count=1;
    FOV_match=[];
    for iSess=1:length(session_nr_vector)
        for iSess2=1:length(session_nr_vector)
            if iSess<iSess2
                im1=ROI_all(iSess).MIP_std;
                im2=ROI_all(iSess2).MIP_std;
                if numel(im1)>numel(im2)
                    A=im1;
                    template=im2;
                else
                    A=im2;
                    template=im1;
                end
                CC=normxcorr2(template,A);
                
                %% get shift coordinates relative to biggest image
                CC_max=max(CC(:));
                [i,j]=find(CC==CC_max);
                
                peakX=j-size(template,2)/2+1;
                peakY=i-size(template,1)/2+1;
                
                FOV_match(count,:)=[count iSess iSess2 CC_max peakX-size(template,2)/2 peakY-size(template,1)/2];
                count=count+1;
                
                if 0
                    %%
                    figure(1)
                    subplot(2,2,1)
                    imshow(im1,[])
                    title(size(im1))
                    subplot(2,2,2)
                    imshow(im2,[])
                    title(size(im2))
                    subplot(2,2,3)
                    imshow(CC,[])
                    subplot(2,2,4)
                    imshow(A,[])
                    title([peakX peakY])
                    hold on
                    plot(peakX,peakY,'m*')
                    plotRect([peakX-size(template,2)/2 peakY-size(template,1)/2 peakX+size(template,2)/2 peakY+size(template,1)/2],'r');
                    hold off
                    colormap(green)
                    
                    pause(.5)
                end
            end
        end
    end
    
    %FOV_match
    
    %% For close distances, we can do nearest neighbor matching
    figure(2)
    count=1;
    min_dist=inf;
    for iSess=1:length(session_nr_vector)
        for iSess2=1:length(session_nr_vector)
            if iSess<iSess2
                
                if 1
                    im1=ROI_all(iSess).MIP_std;
                    im2=ROI_all(iSess2).MIP_std;
                    offset=FOV_match(count,5:6);
                    count=count+1;
                    A=cat(1,ROI_all(iSess).ROI_definitions.center_coords);
                    B=cat(1,ROI_all(iSess2).ROI_definitions.center_coords);
                    
                    B=B+repmat(offset([1 2]),size(B,1),1);
                end
                
                if 0
                    %% Verify match
                    plot(A(:,1),A(:,2),'ro')
                    hold on
                    plot(B(:,1),B(:,2),'bo')
                    hold off
                    title([size(A,1) size(B,1)])
                    axis ij
                    pause(.5)
                end
                
                if 1
                    for iA=1:size(A,1)
                        %%% calc nearest distance between coord
                        ref_coord=A(iA,:)-offset;
                        [min_value,pos]=min(calc_dist(ref_coord,B));
                        
                        match_matrix(iA,:)=[iA pos min_value];
                    end
                    sel=match_matrix(:,3)<min_dist;
                    match_matrix_close=match_matrix(sel,:)
                    
                end
                
                if 1
                    [D,Z,transform]=procrustes(A,B)
                end
                
                if 1
                    %% show matched pairs
                    N=size(match_matrix_close,1);
                    clf
                    hold on
                    for iROI=1:N
                        c1=A(match_matrix_close(iROI,1),:);
                        c2=B(match_matrix_close(iROI,2),:);
                        plot([c1(1) c2(1)],[c1(2) c2(2)],'ro-')
                        plot(c1(1),c1(2),'bo')
                    end
                    hold off
                    pause(.5)
                end
                
                % this matrix represents ROI coords in both FOVs that are
                % so closely matched we can call them the same
            end
        end
    end
    
end






