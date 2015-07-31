classdef imaging_datasets < handle
    properties
        cluster_nr=[];
        session_vector=[];
        nFrames=[];
        FOV_info=struct;
        MIP_std=[];
        ROI_definitions=struct;
        ROI_definition_nr=[];
        STIM=[];
        RESP=[];
        SPIKE=[];
        timeline=[];        
    end
    
    methods
        % Constructor
        function self=imaging_datasets(varargin)
            builder=varargin{1};
            self.ROI_definition_nr=varargin{2};
            self.FOV_info=builder.FOV_info;
            self.MIP_std=builder.MIP_std;
            self.ROI_definitions=builder.ROI_definitions(self.ROI_definition_nr).ROI;            
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
        
        function RF_analysis(varargin)
            self=varargin{1};
            
            if nargin>=2
                iROI=varargin{2};
            else
                iROI=1;
            end
            condition_vector=self.STIM(:,6);
            conditions=unique(condition_vector(condition_vector>-1));
            nConditions=length(conditions);
            
            condition_matrix=parse_conditions(condition_vector);
            nTrials=size(condition_matrix,1); % ignore final trial, because it is rarely complete
            
            frame_selector_trace=[-6 12];
            frame_selector_resp=[2 7];
            calcium_matrix=zeros(nTrials,diff(frame_selector_trace)+1);
            for iTrial=2:nTrials-1
                trial_info=condition_matrix(iTrial,:);
                
                frame_range=trial_info(2)+frame_selector_resp(1):trial_info(2)+frame_selector_resp(2);
                condition_matrix(iTrial,end)=mean(self.RESP(frame_range,iROI));
                
                frame_range=trial_info(2)+frame_selector_trace(1):trial_info(2)+frame_selector_trace(2);
                calcium_matrix(iTrial,:)=self.RESP(frame_range,iROI);
            end
            
            M=pivotTable(condition_matrix,5,'mean',7);
            S=pivotTable(condition_matrix,5,'std',7);
            E=pivotTable(condition_matrix,5,'ste',7);
            %[M S E]
            
            y_range=[-5 40];
            figure(2)
            subplot(211)
            X=self.timeline;
            plot(X,self.RESP(:,iROI))
            axis([X([1 end])' y_range])
            subplot(212)
            if nConditions==32
                MAP=flipud(reshape(M,4,8));
                imagesc(MAP)
                set(gca,'CLim',[-3 3])
                %self.Results.RF_maps(:,:,iROI)=MAP;
            else
                bar(M)
                hold on
                errorbar(M,E,'r.')
                hold off
            end
            
            figure(3)
            nCols=ceil(sqrt(nConditions));
            nRows=ceil(nConditions/nCols);
            for iCondition=1:nConditions
                condition_nr=conditions(iCondition);
                sel=condition_matrix(:,5)==condition_nr;
                avg_trace=mean(calcium_matrix(sel,:));
                subplot(nRows,nCols,iCondition)
                X_AX=frame_selector_trace(1):frame_selector_trace(2);
                cla
                hold on
                y_range=[-1 4];
                plot([0 0],y_range,'r')
                plot([7 7],y_range,'k')
                plot(X_AX,avg_trace)
                
                axis([X_AX([1 end]) y_range])
                title(condition_nr)
                set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            end
            
        end
    end
end