clear all
clc

%% Original file
load('/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-07_AF11/data_analysis/20150407_AF11_006.mat','session_data')
orig=session_data;

%% New file
load('/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-07_AF11_compare/data_analysis/20150407_AF11_006.mat','session_data')
new=session_data;


%% check if all condition info is the same: exclude extraction problem
A=orig.stimulus_matrix_ext;
B=new.Experiment_info.stim_matrix;

%diff([A(:,5) B(:,5)],[],2)
%diff([A(:,8) B(:,6)],[],2)
% nope!

%% check extracted response profiles
A=orig.activity_matrix(:,4);
B=new.Activity_traces.activity_matrix(:,4);

%rescale_factor=max(A)/max(B);
%B=B*rescale_factor;
plot(A)
hold on
plot(B)
hold off
title(corr(A,B))
% traces are really similar, but not identical

%% Check shift_matrix
A=orig.motion_correction.shift_matrix(:,2:3);
B=new.motion_correction.shift_matrix(:,2:3);
[A B]



%% Check average response to positions
condition_vector=new.Experiment_info.stim_matrix(:,6);
condition_matrix_01=parse_conditions(condition_vector);
condition_matrix_02=condition_matrix_01;
conditions=unique(condition_vector(condition_vector>-1));
nConditions=length(conditions);
nTrials=size(condition_matrix_01,1);
frame_selector_trace=[-6 12];
frame_selector_resp=[2 7];
calcium_matrix_01=zeros(nTrials,diff(frame_selector_trace)+1);
calcium_matrix_02=zeros(nTrials,diff(frame_selector_trace)+1);
data_01=orig.activity_matrix(:,4);;
data_02=new.Activity_traces.activity_matrix(:,4);
for iTrial=2:nTrials-1
    trial_info=condition_matrix_01(iTrial,:);
    
    frame_range=trial_info(2)+frame_selector_resp(1):trial_info(2)+frame_selector_resp(2);
    condition_matrix_01(iTrial,end)=mean(data_01(frame_range,end));
    condition_matrix_02(iTrial,end)=mean(data_02(frame_range,end));
    
    frame_range=trial_info(2)+frame_selector_trace(1):trial_info(2)+frame_selector_trace(2);
    calcium_matrix_01(iTrial,:)=data_01(frame_range,end);
    calcium_matrix_02(iTrial,:)=data_02(frame_range,end);
end

figure(3)
nCols=ceil(sqrt(nConditions));
nRows=ceil(nConditions/nCols);
for iCondition=1:nConditions
    condition_nr=conditions(iCondition);
    sel=condition_matrix_01(:,5)==condition_nr;
    avg_trace_01=mean(calcium_matrix_01(sel,:));
    avg_trace_02=mean(calcium_matrix_02(sel,:));
    subplot(nRows,nCols,iCondition)
    X_AX=frame_selector_trace(1):frame_selector_trace(2);
    cla
    hold on
    y_range=[-1 10];
    plot([0 0],y_range,'r')
    plot([7 7],y_range,'k')
    plot(X_AX,avg_trace_01)
    plot(X_AX,avg_trace_02)
    
    axis([X_AX([1 end]) y_range])
    title(condition_nr)
end

figure()
M1=pivotTable(condition_matrix_01,5,'mean',7);
M2=pivotTable(condition_matrix_02,5,'mean',7);
bar([M1 M2])


