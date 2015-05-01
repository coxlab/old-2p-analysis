clear all
clc

header_script

data_folder=['/Users/' user_name '/Dropbox (coxlab)/2p-data/2015-04-15_AF11'];
session_name='20150415_AF11_RM_006.mat';
loadName=fullfile(data_folder,'data_analysis',session_name);

plotIt=3;
plot_individual=0;
nSweeps=2;
fold_over_sweeps=0;

number_of_frames_per_condition=69;
remove_threshold=10;

session_name=strrep(session_name,'_',' ');

if exist(loadName,'file')
    load(loadName,'session_data')
    
    %%
    nROI=length(session_data.ROI_definitions);
    ROI_vector=cat(1,session_data.ROI_definitions.ROI_nr);
    stimulus_matrix_ext=session_data.stimulus_matrix_ext;
    activity_matrix=session_data.activity_matrix;
    
    nFrames=size(activity_matrix,1);
    frameRate=session_data.data(5);
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
    conditions=unique(condition_vector(condition_vector>0));
    nConditions=length(conditions);
    
    transition_points_start=[1; find(diff(condition_vector)>0)+1];
    transition_points_end=[find(diff(condition_vector)<0) ; length(condition_vector)];
    
    %%% clean up
    %[length(transition_points_start) length(transition_points_end)]
    
    %diff(transition_points_start)
    del=find(diff(transition_points_start)<=remove_threshold)+1;
    transition_points_start(del)=[];
    
    %diff(transition_points_end)
    del=find(diff(transition_points_end)<=remove_threshold)+1;
    transition_points_end(del)=[];
    
    %%
    if length(transition_points_start)~=length(transition_points_end)
        [length(transition_points_start) length(transition_points_end)]
        error('Problem parsing out condition start and stop times')
    end
    condition_matrix=[transition_points_start transition_points_end condition_vector(transition_points_start+5) transition_points_end-transition_points_start+1];
    condition_matrix(condition_matrix(:,4)<60,:)=[]; % remove incomplete trials
    nRepeats_per_condition=hist(condition_matrix(:,3),4);
    
    nBlocks_complete=(size(condition_matrix,1)-mod(size(condition_matrix,1),nConditions))/nConditions;
    order_matrix=reshape(condition_matrix(1:nBlocks_complete*nConditions,3),nConditions,[]);
    if any(std(sort(order_matrix),[],2))
        sort(order_matrix)
        error('Condition labels were not read out correctly, go back to step 1a')
    end
    
    condition_matrix_complete=condition_matrix(1:nBlocks_complete*nConditions,:);
    if std(condition_matrix_complete(:,4))~=0
        error('Different number of samples per condition')
    end
    nSamples_per_repeat=mean(condition_matrix_complete(:,4));
    fprintf('Number of repeats per condition: %d\n',nBlocks_complete)
    
    
    %%
    % potential cool cells 3#2 3#6 5#29 6#30(cell)
    [sorted,order]=sort(std(activity_matrix),'descend');
    top_five=order(1:5);
    
    ROI_choice=[];
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
            x1=condition_matrix_complete(iTrial,1);
            x2=condition_matrix_complete(iTrial,2);
            c=condition_matrix_complete(iTrial,3);
            
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
                    cond_info=condition_matrix_complete(condition_matrix_complete(:,3)==iCond,:);
                    cond_matrix=zeros(nBlocks_complete,nSamples_per_repeat);
                    for iBlock=1:nBlocks_complete
                        sel=cond_info(iBlock,1):cond_info(iBlock,2);
                        T=stimulus_matrix_ext(sel,1);
                        T=T-T(1); % fold times
                        trace=activity_matrix(sel,iROI);
                        
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
                            M=mean(cond_matrix);
                            E=ste(cond_matrix);
                            figure(3)
                            subplot(5,1,iCond)
                            plot(T([1 end]),[0 0],'k')
                            hold on
                            %shadedErrorBar(T,M,E,'r',1);
                            shadedErrorBar(T,M,E,colors{iCond},1);
                            hold off
                            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
                            axis([T([1 end])' -5 30])
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
                            axis([T([1 end/nSweeps])' -5 30])
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
end