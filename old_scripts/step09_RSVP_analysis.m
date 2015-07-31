clear all
clc

%%% for the most responsive position, do selectivity analysis over each
%%% cell or each population. Try also to see how generalization works over
%%% all other positions

%%% BV20150518: unbiased method: every
% find number of presentations per unique combination of shape and position
% i.e. 1/384 and concat all the responses per FOV

header_script
%save_it=0;
dataset_selector=3;

nConditions=12;
frame_selector=[2 7];

for data_type=1:2
    sample_interval_seconds=[-2 4];
    switch data_type
        case 1
            y_range=[-5 25];
        case 2
            y_range=[-1 4];
    end
    
    %%% Load requested merged dataset
    loadName=fullfile(data_folder,'data_analysis',sprintf('dataset_%03d.mat',dataset_selector));
    load(loadName,'dataset')
    
    dataset.session_vector'
    sample_rate=1/mean(diff(dataset.timescale));
    
    sample_interval=round(sample_interval_seconds*sample_rate)
    nSamples=range(sample_interval)+1;
    x_range=sample_interval;
    
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
    
    % 0305: 6(11,12,13,14,23) or 7(3,11,24,37) IARPA
    % 0407: 3(4,46,51,53,61,63,65)
    
    switch exp_name
        case '2015-03-05_AF11'
            switch dataset_selector
                case 1
                    ROI_vector=22;            
                case 6
                    ROI_vector=13;%[11 12 13 14 23];
                case 7
                    ROI_vector=[3 11 24 37];
            end
        case '2015-04-07_AF11'
            switch dataset_selector
                case 1
                    ROI_vector=[34];
                case 2
                    ROI_vector=[27 32 44 46 54];
                case 3
                    ROI_vector=4;%[4 46 51 53 61 63 65];
            end
        otherwise
            ROI_vector=1:nROI;
    end
    
    
    for iROI=ROI_vector
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
                    [sorted,order]=sort(vector,'descend');
                    most_active_position=order(3); % => this will lead to an invariance question: how are responses on separate positions related
            end
            
            V=zeros(32,1);
            V(most_active_position)=1;
            RF_disp=reshape(V,4,8)
            
            %%% Get stimulus conditions presented at selected positions
            A=parse_conditions(stim_matrix(:,8)); % position
            B=parse_conditions(stim_matrix(:,5)); % shape
            sel=ismember(A(:,5),most_active_position);
            start_rows=B(sel,2);
            start_rows_pos=A(sel,5);
            stimulus_vector=B(sel,5);
            counts=hist(stimulus_vector,0:11);
            
            response_per_stimulus=zeros(nConditions,1);
            response_per_stimulus_std=response_per_stimulus;
            
            figure(1)
            clf
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
                        elseif repeat_start>nFrames-sample_interval(2)
                            raw_trace=[Y(repeat_start+sample_interval(1):repeat_start);zeros(sample_interval(2),1)];
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
                            %
                            hold on
                            %plot(T,trace_matrix)
                            if nRepeats>1
                                shadedErrorBar(T,mean(trace_matrix,1),ste(trace_matrix,1));
                                CC=corr(trace_matrix');
                                title(sprintf('r=%3.2f (N=%d)',[CC(1,2) nRepeats]))
                            else
                                plot(T,mean(trace_matrix,1))
                            end
                            plot([0 0],y_range,'r-')
                            plot([0 0]+2*sample_rate,y_range,'k-')
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
            
            
            %%% Make figures
            figure(2)
            ROI_definitions=get_ROI_definitions(dataset,ROI_definition_nr);
            ROI_prop=ROI_definitions(iROI);
            coords=cat(1,ROI_prop.center_coords);
            nPix=sum(ROI_prop.mask_soma(:));
            %ROI_prop.ellipse_properties
            [sorted,order]=sort(response_per_stimulus,'descend');
            
            % show cell location on MIP
            subplot(211)
            imshow(calc_gamma(dataset.MIP.data,.3),[])
            hold on
            plot(coords(:,1),coords(:,2),'ro','markersize',30)
            hold off
            colormap(green)
            
            % show receptive field positions used
            subplot(212)
            imshow(imresize(RF_disp,50),[])
            colormap(green)
            
            % show average traces per condition
            figure(3)
            subplot(311)
            plot([start_rows start_rows]',repmat(y_range,length(start_rows),1)','m-')
            hold on
            plot(Y)
            hold off
            box off
            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            
            % show number of repeats per condition
            subplot(312)
            bar(counts(order))
            axis([0 13 0 4])
            box off
            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            set(gca,'xticklabel',order)
            
            % show ranked response per stimulus
            subplot(313)
            bar(response_per_stimulus(order))
            hold on
            errorbar(response_per_stimulus(order),response_per_stimulus_std(order),'r.')
            hold off
            
            axis([0 13 y_range])
            box off
            set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            set(gca,'xticklabel',order)
        end
        
        if save_it
            %%
            %
            save_folder=fullfile(data_folder,'data_analysis','RSVP_selectivity');
            if exist(save_folder,'dir')==0
                mkdir(save_folder)
            end
            cd(save_folder)
            if data_type==1
                saveName=sprintf([exp_name '_%d_%d_traces.eps'],[dataset_selector ROI_prop.ROI_nr]);
            else
                saveName=sprintf([exp_name '_%d_%d_traces_deconv.eps'],[dataset_selector ROI_prop.ROI_nr]);
            end
            print(1,saveName,'-depsc')
            
            %
            saveName=sprintf([exp_name '_%d_%d_ID.eps'],[dataset_selector ROI_prop.ROI_nr]);
            print(2,saveName,'-depsc')
            
            %
            if data_type==1
                saveName=sprintf([exp_name '_%d_%d_overview.eps'],[dataset_selector ROI_prop.ROI_nr]);
            else
                saveName=sprintf([exp_name '_%d_%d_overview_deconv.eps'],[dataset_selector ROI_prop.ROI_nr]);
            end
            print(3,saveName,'-depsc')
        end        
    end
end

