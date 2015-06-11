clear all
clc

header_script

switch 1
    case 1
        data_folder=['/Users/' user_name '/Dropbox (coxlab)/2p-data/2015-04-15_AF11'];
        %session_name='20150415_AF11_RM_006.mat'; % has cool cell
        session_name='20150415_AF11_RM_002.mat'; 
    case 2
        session_name='dataset_001.mat';
end
loadName=fullfile(data_folder,'data_analysis',session_name);


plotIt=3;
plot_individual=0;
nSweeps=2;
fold_over_sweeps=0;

number_of_frames_per_condition=69;
remove_threshold=10;

session_name=strrep(session_name,'_',' ');

if exist(loadName,'file')
    
    load(loadName,'dataset')
    
    if exist('dataset','var')
        disp('Using dataset')
    else
        load(loadName,'session_data')
        disp('Using session_data')
        dataset.stim_matrix=session_data.stimulus_matrix_ext;
        dataset.resp_matrix=session_data.activity_matrix;
        dataset.resp_matrix_NND=session_data.spike_matrix;
        
        dataset.timescale=(1:session_data.data(2))/session_data.data(5)';
        dataset.ROI_definitions=session_data.ROI_definitions;
    end
    
    
    %%
    ROIs=get_ROI_definitions(dataset,ROI_definition_nr);
    nROI=length(ROIs);
    ROI_vector=cat(1,ROIs.ROI_nr);
    stimulus_matrix_ext=dataset.stim_matrix;
    switch 1
        case 1 % regular delta_f/f data
            activity_matrix=dataset.resp_matrix;
        case 2 % use deconvolved data
            activity_matrix=dataset.resp_matrix_NND;
    end
    
    
    nFrames=size(activity_matrix,1);
    %frameRate=session_data.data(5);
    frameRate=1/mode(diff(dataset.timescale));
    time_line=(1:nFrames)/frameRate;
    
    %% look at crosscorrelations
    figure(1)
    clf
    CC=corr(activity_matrix);
    CC_between=CC(eye(nROI)==0);
    [r,c]=find(CC==max(CC_between));
    [ROI_vector(r(1)) ROI_vector(c(1)) max(CC_between)]
    imagesc(CC)
    axis square
    colorbar
    set(gca,'Clim',[-1 1])
    
    %% Check number of repeats per condition
    
    condition_vector=stimulus_matrix_ext(:,5);
    condition_matrix_complete=parse_conditions(condition_vector);
    
    condition_vector=stimulus_matrix_ext(:,5);
    conditions=unique(condition_vector(condition_vector>0));
    nConditions=length(conditions);
    
    %%% Get number of frames for each condition and delete trials
    %%% that have less frames
    %nSamples_per_repeat=mean(condition_matrix_complete(:,4))-1;
    nSamples_per_repeat=mode(condition_matrix_complete(:,4));
    del_incomplete_trials=1;
    if del_incomplete_trials==1
        sel=condition_matrix_complete(:,4)<nSamples_per_repeat;
        condition_matrix_complete(sel,:)=[];
    end
    
    %%% Make nSamples_per_repeat even, since we want to fold the
    %%% two sweeps later.
    if mod(nSamples_per_repeat,2)==1
        nSamples_per_repeat=nSamples_per_repeat-1;
    end
    
    counts=hist(condition_matrix_complete(:,5),1:4)
    nBlocks_complete=min(counts);
    
    %%
    % potential cool cells 3#2 3#6 5#29 6#30(cell)
    [sorted,order]=sort(std(activity_matrix),'descend');
    top_five=order(1:5);
    
    ROI_choice=6;
    if isempty(ROI_choice)
        ROI_choice=top_five(1);
    else
        ROI_choice=find(ROI_vector==ROI_choice);
    end
    %for iROI=
    for iROI=ROI_choice
        ROI_nr=ROI_vector(iROI);
        
        %% Plot raw trace and conditions
        colors={'r','g','m','c'}; % reddish colors for
        y_offset=[-3 -3 -4 -4]; % vertical or horizontal
        figure(2)
        plot(time_line([1 end]),[0 0],'k')
        hold on
        plot(time_line,activity_matrix(:,iROI))
        for iTrial=1:size(condition_matrix_complete,1)
            x1=condition_matrix_complete(iTrial,2);
            x2=condition_matrix_complete(iTrial,3);
            c=condition_matrix_complete(iTrial,5);
            
            plot(time_line([x1 x2]),[y_offset(c) y_offset(c)],colors{c},'lineWidth',3)
        end
        hold off
        axis([time_line([1 end]) -5 20])
        box off
        title(ROI_nr)
        
        %% Analyse per condition
        data_matrix=zeros(nConditions,4);
        for iCond=1:nConditions
            % select all trials for this condition
            condition_nr=conditions(iCond);
            
            sel=stimulus_matrix_ext(:,5)==condition_nr;
            condition_data=stimulus_matrix_ext(sel,:);
            
            data_matrix(iCond,:)=[iCond mean(activity_matrix(sel,iROI)) std(activity_matrix(sel,iROI)) ste(activity_matrix(sel,iROI))];
            
            %%% Show plots
            switch plotIt
                case 0
                case 1
                    %%
                    figure(3)
                    subplot(5,1,iCond)
                    plot(time_line,activity_matrix(:,iROI))
                    hold on
                    %plot(condition_data(:,1),condition_data(:,6:7))
                    plot(stimulus_matrix_ext(sel,1),stimulus_matrix_ext(sel,6)/max(stimulus_matrix_ext(sel,6))-10,'rs')
                    plot(stimulus_matrix_ext(sel,1),stimulus_matrix_ext(sel,7)/max(stimulus_matrix_ext(sel,7))-10,'gs')
                    hold off
                    axis([0 500 -20 20])
                    set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
                case 2
                    X1=stimulus_matrix_ext(sel,6);
                    X2=stimulus_matrix_ext(sel,7);
                    Y=activity_matrix(sel,iROI);
                    
                    figure(3)
                    subplot(5,2,(iCond-1)*2+1)
                    plot(X1,Y,'.')
                    subplot(5,2,(iCond-1)*2+2)
                    [(iCond-1)*2+1 (iCond-1)*2+2]
                    plot(X2,Y,'.')
                    set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
                case 3
                    cond_info=condition_matrix_complete(condition_matrix_complete(:,5)==iCond,:);
                    cond_matrix=zeros(nBlocks_complete,nSamples_per_repeat);
                    for iBlock=1:nBlocks_complete
                        sel=cond_info(iBlock,2):cond_info(iBlock,3);
                        T=stimulus_matrix_ext(sel,1);
                        T=T(1:nSamples_per_repeat)-T(1); % fold times
                        trace=activity_matrix(sel,iROI);
                        trace=trace(1:nSamples_per_repeat);
                        
                        cond_matrix(iBlock,:)=trace;
                        
                        cond_matrix_folded=reshape(cond_matrix',size(cond_matrix,2)/nSweeps,[])';
                        if plot_individual==1
                            %%
                            figure(3)
                            subplot(5,1,iCond)
                            hold on
                            plot(T,trace)
                            hold off
                        end
                    end
                    
                    if plot_individual==0
                        if fold_over_sweeps==0
                            M=mean(cond_matrix,1)';
                            E=ste(cond_matrix,1)';
                            figure(3)
                            subplot(5,1,iCond)
                            plot(T([1 end]),[0 0],'k')
                            hold on
                            %shadedErrorBar(T,M,E,'r',1);
                            shadedErrorBar(T(:),M(:),E(:),colors{iCond},1);
                            hold off
                            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
                            %axis([T([1 end])' -5 30])
                        else
                            M=mean(cond_matrix_folded);
                            E=ste(cond_matrix_folded);
                            figure(3)
                            subplot(5,1,iCond)
                            plot(T([1 end/nSweeps]),[0 0],'k')
                            hold on
                            %shadedErrorBar(T,M,E,'r',1);
                            shadedErrorBar(T(1:end/nSweeps),M,E,colors{iCond},1);
                            hold off
                            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
                            %axis([T([1 end/nSweeps])' -5 30])
                        end
                    end
            end
        end
        
        % Show average traces
        subplot(515)
        bar(data_matrix(:,2))
        hold on
        errorbar(data_matrix(:,2),data_matrix(:,4),'r.')
        hold off
        axis([0 5 -.5 7])
        box off
        title(session_name)
        xlabel(ROI_nr)
        set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
    end
    
    
    %% Analyse all ROIs
    ROI_data=struct;
    for iROI=1:nROI
        data_matrix=zeros(nConditions,4);
        for iCond=1:nConditions
            % select all trials for this condition            
            condition_nr=conditions(iCond);
            sel=condition_matrix_complete(:,5)==condition_nr;
            cond_matrix=condition_matrix_complete(sel,:);
            nRepeats=size(cond_matrix,1);
            
            traces=[];
            for iRepeat=1:nRepeats
                idx=cond_matrix(iRepeat,2):cond_matrix(iRepeat,3);
                trace=activity_matrix(idx,iROI);
                traces=cat(2,traces,trace);
            end
            ROI_data(iROI).conditions(iCond).traces=traces;
            ROI_data(iROI).conditions(iCond).avg_trace=mean(traces,2);
            ROI_data(iROI).conditions(iCond).avg_trace_smooth=mean(traces,2);
            ROI_data(iROI).conditions(iCond).std_trace=std(traces,[],2);
            ROI_data(iROI).conditions(iCond).ste_trace=ste(traces,2);
            % add stuff like, where is max? and FFT analysis            
            
            
        end
        plot_it=1;
        if plot_it==1
            nCols=ceil(sqrt(nROI));
            nRows=ceil(nROI/nCols);
            trace_all=cat(2,ROI_data(iROI).conditions.avg_trace);
            subplot(nRows,nCols,iROI)
            plot(trace_all)
            axis([0 60 -5 30])
            title(ROIs(iROI).ROI_nr)
            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
        end
    end
    
end