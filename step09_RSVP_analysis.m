clear all
clc

%%% for the most responsive position, do selectivity analysis over each
%%% cell or each population. Try also to see how generalization works over
%%% all other positions

%%% BV20150518: unbiased method: every
% find number of presentations per unique combination of shape and position
% i.e. 1/384 and concat all the responses per FOV

header_script
save_it=0;
dataset_selector=2; % 6(13,23) or 7(11,37) IARPA
data_type=1;

nConditions=12;
frame_selector=[2 7];

%%% Load requested merged dataset
loadName=fullfile(data_folder,'data_analysis',sprintf('dataset_%03d.mat',dataset_selector));
load(loadName,'dataset')

dataset.session_vector'

nFrames=dataset.nFrames;
nROI=dataset.nROI;
stim_matrix=dataset.stim_matrix;
resp_matrix=dataset.resp_matrix;
resp_matrix_NND=dataset.resp_matrix_NND;

%% From the RF MAPs, find significant positions
nCols=ceil(sqrt(nROI));
nRows=ceil(nROI/nCols);
TH=0.8;
most_active_positions=NaN(nROI,1);
for iROI=1%1:nROI
    RF=flipud(dataset.MAPs(iROI).MAP_zscored);
    vector=RF(:);
    sel=vector>TH;
    
    switch data_type
        case 1
            Y=resp_matrix(:,iROI);
        case 2
            Y=resp_matrix_NND(:,iROI);
    end    
    
    if sum(sel)==0
        figure(1)
        clf
    else
        switch 1
            case 0
                most_active_position=1:32;
            case 1 % safe
                nPositions=min([sum(sel) 3]);
                [sorted,order]=sort(vector,'descend');
                most_active_position=order(1:nPositions);
            case 2 % unsafe
                nPositions=3;
                [sorted,order]=sort(vector,'descend');
                most_active_position=order(1:nPositions);
            case 3
                nPositions=8;
                [sorted,order]=sort(vector,'descend');
                most_active_position=order(1:nPositions);
            case 4
                most_active_position=[22 26];
            case 5
                nPositions=min([sum(sel) 1]);
                [sorted,order]=sort(vector,'descend');
                most_active_position=order(1);
        end
        
        V=zeros(32,1);
        V(most_active_position)=1;
        reshape(V,4,8)
        
        sample_interval=[-20 60];
        nSamples=range(sample_interval)+1;
        x_range=sample_interval;
        y_range=[-5 20];
        
        %%% Get stimulus conditions presented at selected positions
        A=parse_conditions(stim_matrix(:,8)); % position
        B=parse_conditions(stim_matrix(:,5)); % shape
        sel=ismember(A(:,5),most_active_position);
        start_rows=B(sel,2);
        stimulus_vector=B(sel,5);
        counts=hist(stimulus_vector,0:11);
        
        response_per_stimulus=zeros(nConditions,1);
        response_per_stimulus_std=response_per_stimulus;
        for iCond=1:nConditions
            stim_nr=iCond-1;
            repeat_starts=start_rows(stimulus_vector==stim_nr);
            nRepeats=counts(iCond);
            if nRepeats>0
                % get average response for each repeat                
                repeat_vector=zeros(nRepeats,1);
                trace_matrix=zeros(nRepeats,nSamples);
                for iRepeat=1:nRepeats
                    repeat_start=repeat_starts(iRepeat);
                    if repeat_start<abs(sample_interval(1))
                        raw_trace=[zeros(abs(sample_interval(1)),1) ; Y(repeat_start:repeat_start+sample_interval(2))];
                    else 
                        raw_trace=Y(repeat_start+sample_interval(1):repeat_start+sample_interval(2));
                    end
                    trace_matrix(iRepeat,:)=raw_trace;
                    repeat_vector(iRepeat)=mean(Y(repeat_start+frame_selector(1):repeat_start+frame_selector(2)));
                    %stim_matrix(repeat_start-2:repeat_start+6,:); % sanity check
                end
                
                figure(1)
                subplot(4,3,iCond)                
                switch 1
                    case 1
                        T=sample_interval(1):sample_interval(2);
                        plot(T,mean(trace_matrix,1))
                        hold on
                        plot([0 0],y_range,'r-')
                        plot([0 0]+24,y_range,'k-')
                        hold off
                        axis([x_range y_range])
                    case 2
                        plot(Y)
                        hold on
                        plot(repeat_starts,ones(size(repeat_starts)),'m*')
                        hold off
                end
                
                response_per_stimulus(iCond)=mean(repeat_vector);
                response_per_stimulus_std(iCond)=ste(repeat_vector);
            else
                disp(['no data for condition ' iCond])
                subplot(4,3,iCond)
                cla
            end
        end
        
        %%
        figure(2)
        subplot(311)
        plot(Y)
        hold on
        plot(start_rows,ones(size(start_rows)),'m*')
        hold off
        subplot(312)
        bar(counts)
        subplot(313)
        [sorted,order]=sort(response_per_stimulus,'descend');
        bar(response_per_stimulus(order))
        hold on
        errorbar(response_per_stimulus(order),response_per_stimulus_std(order),'r.')
        hold off
        
    end
    
    if plot_it==1
        subplot(nRows,nCols,iROI)
        imshow(RF>TH)
        title(iROI)
    end
    
    
end


%
% %%
% stim_matrix_crop=stim_matrix;
% stim_matrix_crop(stim_matrix_crop(:,4)==-1,:)=[];
% M=[pivotTable(stim_matrix_crop,8,'mean',8) pivotTable(stim_matrix_crop,8,'mean',5) pivotTable(stim_matrix_crop,8,'length',8)];
%
%
% M
% size(M)
%
% %% Collect average responses for all ROIs
%
%
%
% %%
%



