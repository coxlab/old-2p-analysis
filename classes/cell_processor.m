%%% Loop over all cells of a certain session per FOV
%%% Copy demographic info from dataset
%%% Add methods to run a certain type of analysis routine and save those
%%% data into the properties
%%% Should be easy to combine certain property over multiple cells


classdef cell_processor < handle
    properties
        %%% Information about animal, day, calibration, session, protocol
        %%% ...
        animal_ID='';
        session_date=[];
        exp_type=[];
        exp_name='';
        FOV_info=struct;
        
        cell_id=[];
        
        trace=[];
        nTrials=[];
        condition_matrix=[];
        response_matrix=[];
        trials=struct;
        nPerm=100;
        stimulus_col_nr=3;
        position_col_nr=4;
        response_col_nr=5;
        condition_info=struct;
        RF_map=[];
        MU=[];
        SIGMA=[];
        RF_map_norm=[];
        
        TH=[];
        RF_map_TH=[];
        
        %RF_results=struct;
        stim_ID_all=[];
        stim_response_all=[];
        sparseness_all=[];
        
        nResponsive_positions=[];
        stim_ID_resp=[];
        stim_response_resp=[];
        sparseness_avg=[];
        
        stim_response_per_position=struct('condition_nr',[],'stim_ID_vector',[],'response',[])
        sparseness_per_position=[];
        invariance_avg=[];
    end
    
    methods
        function self=cell_processor(varargin)
            %%% Give
            self.cell_id=varargin{1};
            dataset=varargin{2};
            self.animal_ID=dataset.animal_ID;
            self.session_date=dataset.session_date;
            self.exp_name=dataset.exp_name;
            self.exp_type=dataset.exp_type;
        end
        
        function build_condition_matrix(varargin)
            self=varargin{1};
            STIM=varargin{2};
            self.condition_matrix=parse_conditions(STIM(:,4));
            A=parse_conditions(STIM(:,5));
            self.condition_matrix(:,6)=A(:,5);
            B=parse_conditions(STIM(:,6));
            self.condition_matrix(:,7)=B(:,5);
            self.nTrials=size(self.condition_matrix,1);
            %hist(condition_matrix(:,5))
        end
        
        function add_trace(varargin)
            self=varargin{1};
            self.trace=varargin{2};
        end
        
        %%% Receptive field mapping
        function do_RF_analysis(varargin)
            self=varargin{1};
            %parameter=varargin{2};
            
            [self.response_matrix,self.trials]=self.calc_response_matrix();
            self.RF_map=reshape(pivotTable(self.response_matrix,self.position_col_nr,'mean',self.response_col_nr),4,8);
        end
        
        function do_RF_analysis_shuffled(varargin)
            self=varargin{1};
            %parameter=varargin{2};
            trace_shuffled=self.trace(randperm(length(self.trace)));
            mu_vector=zeros(self.nPerm,1);
            sigma_vector=mu_vector;
            for iPerm=1:self.nPerm
                RM_shuffled=self.calc_response_matrix(trace_shuffled);
                RF_map_shuffled=reshape(pivotTable(RM_shuffled,self.position_col_nr,'mean',self.response_col_nr),4,8);
                mu_vector(iPerm)=mean(RF_map_shuffled(:));
                sigma_vector(iPerm)=std(RF_map_shuffled(:));
            end
            
            % Generate normalized map (could be separate function, not for now)
            self.MU=mean(mu_vector);
            self.SIGMA=mean(sigma_vector);
            self.RF_map_norm=(self.RF_map-self.MU)/self.SIGMA;
        end
        
        function [RM,TR]=calc_response_matrix(varargin)
            self=varargin{1};
            % if omitted, using actual trace, if not, allows for running the same analysis on modified trace
            if nargin>=2&&~isempty(varargin{2})
                TRACE=varargin{2};
            else
                TRACE=self.trace;
            end
            
            TR=struct; % clear trials
            for iTrial=1:self.nTrials
                trial_data=self.condition_matrix(iTrial,:);
                %trace(trial_data(2):trial_data(3))
                TR(iTrial).trial_nr=iTrial;
                TR(iTrial).condition_nr=trial_data(5);
                TR(iTrial).stim_nr=trial_data(6);
                TR(iTrial).pos_nr=trial_data(7);
                TR(iTrial).frames=trial_data(2)+1:trial_data(3)+3;
                TR(iTrial).nFrames=length(TR(iTrial).frames);
                TR(iTrial).response=TRACE(TR(iTrial).frames);
                TR(iTrial).response_avg=mean(TR(iTrial).response);
                TR(iTrial).response_std=std(TR(iTrial).response);
            end
            RM=[cat(1,TR.trial_nr) cat(1,TR.condition_nr) cat(1,TR.stim_nr) cat(1,TR.pos_nr) cat(1,TR.response_avg) cat(1,TR.response_std)];
        end
        
        function do_threshold(varargin)
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                self.TH=varargin{2};
            else
                self.TH=1;
            end
            
            self.RF_map_TH=self.RF_map>self.TH;
        end
        
        
        function show_RF_map(varargin)
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                MAP=varargin{2};
            else
                MAP=self.RF_map;
            end
            
            if std(MAP)>eps
                imagesc(MAP)
                max_val=max(abs(MAP(:)));
                range=[-1 1]*max_val;
                set(gca,'cLim',range)
                axis xy
                colorbar
            end
        end
        
        %%% Stimulus selectivity
        function do_stimSelect_analysis(varargin)
            self=varargin{1};
            %parameter=varargin{2};
            
            [RM,TR]=self.calc_response_matrix();
            
            %% over all positions
            self.stim_ID_all=pivotTable(RM,self.stimulus_col_nr,'mean',self.stimulus_col_nr);
            self.stim_response_all=pivotTable(RM,self.stimulus_col_nr,'mean',self.response_col_nr);
            
            %% find significant positions
            position_nr_vector=find(self.RF_map_TH(:)==1);
            if ~isempty(position_nr_vector)
                %% over all responsive positions
                RM_responsive_positions=RM(ismember(RM(:,self.position_col_nr),position_nr_vector),:);
                self.stim_ID_resp=pivotTable(RM_responsive_positions,self.stimulus_col_nr,'mean',self.stimulus_col_nr);
                self.stim_response_resp=pivotTable(RM_responsive_positions,self.stimulus_col_nr,'mean',self.response_col_nr);
                
                self.nResponsive_positions=length(position_nr_vector);
                for iPos=1:self.nResponsive_positions
                    condition_nr=position_nr_vector(iPos);
                    RM_pos=RM(RM(:,self.position_col_nr)==condition_nr,:);
                    stim_ID_vector=pivotTable(RM_pos,self.stimulus_col_nr,'mean',self.stimulus_col_nr);
                    stim_resp=pivotTable(RM_pos,self.stimulus_col_nr,'mean',self.response_col_nr);
                    
                    self.stim_response_per_position(iPos).condition_nr=condition_nr;
                    self.stim_response_per_position(iPos).stim_ID_vector=stim_ID_vector;
                    self.stim_response_per_position(iPos).response=stim_resp;
                end
            end
        end
        
        function get_sparseness(varargin)
            self=varargin{1};
            
            self.sparseness_all=self.calc_sparseness(self.stim_response_all);
            self.sparseness_avg=self.calc_sparseness(self.stim_response_resp);
            
            if self.nResponsive_positions>0
                for iPos=1:self.nResponsive_positions
                    self.sparseness_per_position(iPos)=self.calc_sparseness(self.stim_response_per_position(iPos).response);
                end
            end
        end
        
        function S=calc_sparseness(varargin)
            % Rust and DiCarlo 2012 => Vinje and Gallant 2002
            A=varargin{2};
            N=length(A);
            A(A<0)=0; % rectify to avoid out of bound values
            a=( (sum(A)/N)^2 ) / ( sum((A.^2)/N ) );
            S=(1-a)/(1-1/N);
        end
        
        function calc_invariance(varargin)
            % Rust and DiCarlo 2012
            self=varargin{1};
            
            if self.nResponsive_positions>1
                pairs=possibleComparisons(self.nResponsive_positions);
                nPairs=size(pairs,1);
                
                angular_similarity_vector=zeros(nPairs,1);
                for iPair=1:nPairs
                    pair=pairs(iPair,:);
                    
                    %%% find stimuli that were presented at least once at each of
                    %%% the two positions
                    X1=self.stim_response_per_position(pair(1)).stim_ID_vector;
                    R1=self.stim_response_per_position(pair(1)).response;
                    
                    X2=self.stim_response_per_position(pair(2)).stim_ID_vector;
                    R2=self.stim_response_per_position(pair(2)).response;
                    
                    common=intersect(X1,X2);
                    A_common=R1(ismember(X1,common));
                    B_common=R2(ismember(X2,common));
                    angular_similarity=1 - acos( sum(A_common.*B_common) / (norm(A_common)*norm(B_common)) )/pi*2;
                    if isnan(angular_similarity)
                        angular_similarity=0;
                    end
                    angular_similarity_vector(iPair)=angular_similarity;
                    %r=corr([R1(ismember(X1,common)) R2(ismember(X2,common))]);
                end
                self.invariance_avg=mean(angular_similarity_vector);
            end
        end
        
        
        %         function g(varargin)
        %             self=varargin{1};
        %         end
    end
end