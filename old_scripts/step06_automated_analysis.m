clear all
clc

% load ROIs from different files and check their alignment.
% additional info to be used
% - FOV position
% - coregistration of MIP in both, can be different size

header_script

folder_selector=1;

switch folder_selector
    case 1
        data_folder=fullfile(data_root,'2015-03-05_AF11/data_analysis/');
        
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
        data_folder=fullfile(data_root,'2015-03-05_AF11/data_analysis/');
        loadName_format='20150304_AF11_%03d.mat';
        cluster_nr=1;
        session_nr_vector=[2 3 5];
end

nSessions=length(session_nr_vector);


%ROI_all=struct;
for iSess=1:nSessions
    session_nr=session_nr_vector(iSess);
    loadName=fullfile(data_folder,sprintf(loadName_format,session_nr));
    
    if exist(loadName,'file')
        load(loadName,'session_data')
        ROI_all(iSess)=session_data;
    end
end


%%
%%% given we have copied ROIs, we can be sure that same numbers indicate
%%% same cells. Deleting ROI and creating new ones results in new unique
%%% numbers.

% get all unique numbers
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
data_matrix_all=[];
for iSess=1:nSessions
    cell_numbers=cat(1,ROI_all(iSess).ROI_definitions.ROI_nr);
    cell_selection=ismember(cell_numbers,unique_cell_numbers);
    time_selection=2:size(ROI_all(iSess).stimulus_matrix_ext,1);
    
    STIM=ROI_all(iSess).stimulus_matrix_ext(time_selection,:);
    RESP=ROI_all(iSess).activity_matrix(time_selection,cell_selection);
    %[size(STIM) size(RESP)]
    STIM_ALL=cat(1,STIM_ALL,STIM);
    RESP_ALL=cat(1,RESP_ALL,RESP);
    
    % chop up RESP into trial responses using STIM
    trial_matrix=STIM(STIM(:,4)>0,[1 4]);
    
    % find frames where new trial begins
    new_trial_frames=[2; trial_matrix([0;diff(trial_matrix(:,2))]==1,1)]; % new trial at start of exp
    
    nTrials=max(trial_matrix(:,2));
    
    data_matrix=zeros(nTrials-2,10+nROI);
    data_matrix_base=data_matrix;
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
        trial_data=RESP(sel,:); % all frames collected during a given trial
        trial_data_base=RESP(sel_base,:); % all frames collected during a given trial
        data_matrix(iTrial-1,:)=[iSess iTrial stim_properties mean(trial_data,1)];
        data_matrix_base(iTrial-1,:)=[iSess iTrial stim_properties mean(trial_data_base,1)];
    end
    
    data_matrix_all=cat(1,data_matrix_all,data_matrix);
    data_matrix_base_all=cat(1,data_matrix_all,data_matrix_base);
end

%size(STIM_ALL)
%size(RESP_ALL)
nFrames=size(RESP_ALL,1);
%%
data_matrix=zeros(nROI,8);
for iROI=1:nROI
    cell_nr=unique_cell_numbers(iROI);
    F=RESP_ALL(:,iROI);
    
    TH=std(F)*6;
    F_TH=F>TH;
    rising_edges=find(diff(F_TH)==1); % find all rising edges in thresholded trace
    
    nSpikes=length(rising_edges);
    
    data_matrix(iROI,1:4)=[iROI cell_nr TH nSpikes];
    
    if 1
        %%
        plot(F)
        hold on
        line([0 nFrames],[TH TH],'color','r')
        plot(rising_edges,ones(nSpikes,1)*TH,'m*')
        hold off
        axis([0 nFrames -1 4])
        box off
        title(nSpikes)
    end
end


nPositions=32;
resp_data=zeros(nPositions,nROI,3);
base_data=resp_data;
Y=[];groups=[];
for iPos=1:nPositions
    sel=data_matrix_all(:,9)==iPos;
    resp_matrix=data_matrix_all(sel,11:end);
    
    sel=data_matrix_base_all(:,9)==iPos;
    base_matrix=data_matrix_base_all(sel,11:end);
    
    %nTrials=size(resp_matrix,1);
    resp_data(iPos,:,:)=cat(3,mean(resp_matrix,1),std(resp_matrix,[],1),ste(resp_matrix,1));
    base_data(iPos,:,:)=cat(3,mean(base_matrix,1),std(base_matrix,[],1),ste(base_matrix,1));
    net_data=cat(3,mean(resp_matrix,1)-mean(base_matrix,1),std(resp_matrix,[],1),ste(resp_matrix,1));
    
    % build anova matrix
    Y=cat(1,Y,resp_matrix);
    groups=cat(1,groups,repmat(iPos,size(resp_matrix,1),1));
end

for iROI=1:nROI
    [p_pos,t,stats]=anovan(Y(:,iROI),groups,'display','off');
    data_matrix(iROI,5)=p_pos;
end

data_matrix=sortrows(data_matrix,-4);


%% cell 25 55 of folder 2 seems to be position selective
iROI=11;
ROI_nr=data_matrix(iROI,2)'
nSpikes=data_matrix(iROI,4)

figure(2)
subplot(221)
bar(resp_data(:,iROI,1))
hold on
%errorbar(resp_data(:,iROI,1),resp_data(:,iROI,2),'r.')
errorbar(resp_data(:,iROI,1),resp_data(:,iROI,3),'r.')
hold off
set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
title('Screen Position')
xlabel(sprintf('Anova p=%3.4f',data_matrix(iROI,5)))
axis([1 nPositions -0.5 max(resp_data(:,iROI,1))*2])

subplot(222)
RF=flipud(reshape(resp_data(:,iROI,1),4,8));
max_val=max(abs([min(RF(:)) max(RF(:))]));
imagesc(RF)
axis equal
axis tight

set(gca,'CLim',[-max_val max_val],'xTickLabel',[],'yTickLabel',[])
title('''RF''')
set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
colormap(parula)

% Examine responses
stim_vector=STIM_ALL(:,4)>0;
resp_vector=RESP_ALL(:,iROI);
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

subplot(2,2,[3 4])
max_val=max([4 max(resp_vector)*1.2]);
bar(T,spherify(stim_vector,2)*max_val,'barWidth',1,'FaceColor',[1 1 1]*.9,'EdgeColor',[1 1 1]*.9)
hold on
bar(T,spherify(best_conditions,2)*-1,'barWidth',1,'FaceColor',[1 0 0]*.5,'EdgeColor',[1 0 0]*.5)
bar(T,spherify(best_condition,2)*-1,'barWidth',1,'FaceColor',[1 0 0]*.9,'EdgeColor',[1 0 0]*.9)
plot(T,resp_vector,'color','k')
hold off
axis([0 total_duration -max_val*.25 max_val])
set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
title(sprintf('Cell #%d',ROI_nr))
xlabel('Time (seconds)')
ylabel('\DeltaF/F')

%%% double check
% sel=data_matrix_all(:,8)==order(2);
% resp_matrix=data_matrix_all(sel,10:end);
% resp_matrix(:,iROI)
% nTrials=size(resp_matrix,1);
% new=cat(3,mean(resp_matrix,1)-1,std(resp_matrix,[],1),ste(resp_matrix,1));
