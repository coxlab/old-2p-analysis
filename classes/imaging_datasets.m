classdef imaging_datasets < handle
    properties
        animal_ID='';
        session_date=[];
        session_nr=[];
        cluster_nr=[];
        session_vector=[];
        nSessions=[];
        nFrames=[];
        frame_rate=[];
        nROIs=[];
        FOV_info=struct;
        MIP_std=[];
        ROI_definitions=struct;
        ROI_definition_nr=[];
        exp_type=[];
        exp_name='';
        stim_duration=[];        
        STIM=[];
        RESP=[];
        SPIKE=[];
        timeline=[];
    end
    
    methods
        % Constructor
        function self=imaging_datasets(varargin)
            if nargin>=1
                builder=varargin{1};
                self.ROI_definition_nr=builder.ROI_definition_nr;
                self.FOV_info=builder.FOV_info;
                self.MIP_std=builder.MIP_std;
                self.ROI_definitions=builder.ROI_definitions(self.ROI_definition_nr).ROI;
                self.nROIs=length(self.ROI_definitions);
            end
        end
        
        % Data analysis functions
        function plot_traces(varargin)
            self=varargin{1};
            A=self.RESP;
            nROI=size(A,2);
            nCols=ceil(sqrt(nROI));
            nRows=ceil(nROI/nCols);
            figure()
            for iROI=1:nROI
                subplot(nRows,nCols,iROI)
                plot(A(:,iROI))
                axis([1 size(A,1) -10 40])
                title(sprintf('ROI #%d',iROI))
                set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            end
        end
        
        function plot_FOV(varargin)
            self=varargin{1};
            
            N=length(self);
            figure(23)
            clf
            hold on
            for iFOV=1:N
                dataset_nr=self(iFOV).cluster_nr;
                session_data=self(iFOV);
                circle([0 0],2,100,'r-',2);
                center=session_data.FOV_info.center;
                FOV_rect=[0 0 session_data.FOV_info.size_um];
                ROI=CenterRectOnPoint(FOV_rect,center(1),center(2))/1000;
                
                plotRect(ROI,'k');
                
                axis square
                axis equal
                
                %title(session_data.exp_name)
                text(center(1)/1000,center(2)/1000,sprintf('Depth %3.1fµm',session_data.FOV_info.Z_depth))
                text(ROI(1),ROI(2)+.1,sprintf('#%d',dataset_nr))
            end
            hold off            
        end
        
        function results=RF_analysis(varargin)
            self=varargin{1};
            
            if nargin>=2&&~isempty(varargin{2})
                ROI_vector=varargin{2};
            else
                ROI_vector=1:self.nROIs;
            end
            
            if nargin>=3&&~isempty(varargin{3})
                deconvolve=varargin{3};
            else
                deconvolve=0;
            end
            
            results=struct;
            
            condition_vector=self.STIM(:,6);
            conditions=unique(condition_vector(condition_vector>-1));
            nConditions=length(conditions);
            
            condition_matrix=parse_conditions(condition_vector);
            nTrials=size(condition_matrix,1); % ignore final trial, because it is rarely complete
            
            frame_selector_trace=[-6 12];
            %frame_selector_resp=[1 4];
            frame_selector_resp=round([0 self.stim_duration(1)/1000]*self.frame_rate);            
            
            results.nROIs=length(ROI_vector);
            results.ROI_vector=ROI_vector;
            results.deconvolve=deconvolve;
            results.condition_matrix=condition_matrix;
            results.condition_vector=conditions;
            results.nConditions=nConditions;
            results.nTrials=nTrials;
            results.frame_selector_trace=frame_selector_trace;
            results.frame_selector_resp=frame_selector_resp;
            
            for iROI=1:length(ROI_vector)
                ROI_nr=ROI_vector(iROI);
                calcium_matrix=zeros(nTrials,diff(frame_selector_trace)+1);
                for iTrial=2:nTrials-1
                    trial_info=condition_matrix(iTrial,:);
                    
                    frame_range=trial_info(2)+frame_selector_resp(1):trial_info(2)+frame_selector_resp(2);
                    if deconvolve==0
                        condition_matrix(iTrial,end)=mean(self.RESP(frame_range,ROI_nr));
                    else
                        condition_matrix(iTrial,end)=mean(self.SPIKE(frame_range,ROI_nr));
                    end
                    
                    frame_range=trial_info(2)+frame_selector_trace(1):trial_info(2)+frame_selector_trace(2);
                    nFrames_sel=length(frame_range);
                    if deconvolve==0
                        calcium_matrix(iTrial,:)=self.RESP(frame_range,ROI_nr);
                    else
                        calcium_matrix(iTrial,:)=self.SPIKE(frame_range,ROI_nr);
                    end
                end
                
                cond_sel=unique(condition_matrix(:,5));
                if length(cond_sel)==32
                    M=pivotTable(condition_matrix,5,'mean',7);
                    S=pivotTable(condition_matrix,5,'std',7);
                    E=pivotTable(condition_matrix,5,'ste',7);
                else
                    M(cond_sel)=pivotTable(condition_matrix,5,'mean',7);
                    S(cond_sel)=pivotTable(condition_matrix,5,'std',7);
                    E(cond_sel)=pivotTable(condition_matrix,5,'ste',7);
                end
                X=self.timeline;
                MAP=(reshape(M,4,8));
                
                results.condAverage(iROI).cond_X=X;
                results.condAverage(iROI).cond_mean=M;
                results.condAverage(iROI).cond_std=S;
                results.condAverage(iROI).cond_err=E;
                
                results.condAverage(iROI).RF_map=MAP;
                
                %frame_rate=1/mean(diff(self.timeline));
                X=(frame_selector_trace(1):frame_selector_trace(2))/self.frame_rate;
                results.condTraces(iROI).nRepeats=zeros(nConditions,1);
                results.condTraces(iROI).avg_X=X;
                results.condTraces(iROI).avg=zeros(nFrames_sel,nConditions);
                results.condTraces(iROI).std=zeros(nFrames_sel,nConditions);
                
                for iCondition=1:nConditions
                    condition_nr=conditions(iCondition);
                    sel=condition_matrix(:,5)==condition_nr;
                    nRepeats=sum(sel);
                    avg_trace=mean(calcium_matrix(sel,:));
                    std_trace=std(calcium_matrix(sel,:));
                    ste_trace=ste(calcium_matrix(sel,:));
                    
                    results.condTraces(iROI).nRepeats(iCondition)=nRepeats;
                    results.condTraces(iROI).trace(iCondition).avg=avg_trace;
                    results.condTraces(iROI).trace(iCondition).std=std_trace;
                    results.condTraces(iROI).trace(iCondition).ste=ste_trace;
                end
            end
        end
        
        function plot_RF_map(varargin)
            self=varargin{1};
            if length(self)>1
                error('Please select a single dataset...')
            end
            
            if nargin>=2&&~isempty(varargin{2})
                ROI_vector=varargin{2};
            else
                ROI_vector=1:self.nROIs;
            end
            
            if nargin>=3&&~isempty(varargin{3})
                deconvolve=varargin{3};
            else
                deconvolve=0;
            end
            
            results=self.RF_analysis(ROI_vector,deconvolve);
            N=results.nROIs;
            
            if N==1
                iROI=ROI_vector(1);
                onset_times=find(diff(self.STIM(:,3))==1);
                offset_times=find(diff(self.STIM(:,3))==-1);

                %%% plot trace
                if results.deconvolve==0
                    y_range=[-5 40];
                else
                    y_range=[-1 4];
                end
                figure(2)
                subplot(211)
                X=self.timeline;
                %bar(X,self.STIM(:,3)*y_range(2))
                plot(X([onset_times onset_times])',y_range,'r')
                hold on
                plot(X([offset_times offset_times])',y_range,'k')
                if results.deconvolve==0
                    plot(X,self.RESP(:,iROI))
                else
                    plot(X,self.SPIKE(:,iROI))
                end
                hold off
                axis([X([1 end])' y_range])
                %%% plot RF map
                subplot(212)
                MAP=results.condAverage.RF_map;
                imagesc(MAP)
                axis xy
                if results.deconvolve==0
                    set(gca,'CLim',[-3 3])
                else
                    set(gca,'CLim',[-1 1]*.5)
                end
                
                %%% plot avg trace per condition
                figure(3)
                clf
                nConditions=32;
                nCols=8;
                nRows=ceil(nConditions/nCols);
                X=results.condTraces.avg_X;
                nRepeats=results.condTraces.nRepeats;
                for iCondition=1:nConditions
                    condition_nr=results.condition_vector(iCondition);
                    %stim_duration=1;
                    
                    subplot(nRows,nCols,condition_nr)
                    cla
                    hold on
                    if deconvolve==0
                        y_range=[-1 10];
                    else
                        y_range=[-1 4];
                    end
                    plot([0 0],y_range,'r')
                    plot([self.stim_duration(1) self.stim_duration(1)]/1e3,y_range,'k')
                    shadedErrorBar(X,results.condTraces.trace(condition_nr).avg,results.condTraces.trace(condition_nr).ste);
                    
                    axis([X([1 end]) y_range])
                    title(sprintf('Cond #%d (N=%d)',[condition_nr nRepeats(condition_nr)]))
                    set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
                end
                
                
            else % plot overview
                figure()
                clf
                nCols=ceil(sqrt(N));
                nRows=ceil(N/nCols);
                for iROI=1:N
                    ROI_nr=results.ROI_vector(iROI);
                    subplot(nRows,nCols,iROI)
                    MAP=results.condAverage(iROI).RF_map;
                    imagesc(MAP)
                    axis xy
                    if results.deconvolve==0
                        set(gca,'CLim',[-3 3])
                    else
                        set(gca,'CLim',[-1 1]*.25)
                    end
                    title(sprintf('#%d',ROI_nr))
                end
            end
        end
    end
end