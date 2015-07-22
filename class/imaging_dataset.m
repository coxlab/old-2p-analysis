classdef imaging_dataset < handle
    properties
        file_name='';
        save_name='';
        folder_info=struct('main_folder','2p-data','root_folder','','rel_path','','data_folder','','raw_name','','save_folder','','save_folder_root','data_analysis')
        mov_info=struct('nFrames',[],'Width',[],'Height',[],'state',[],'frame_rate',[],'mov_start_time',[],'mean_lum',[],'blank_frames',[])
        frame_info=struct('bitCode_vector',[],'main_bitCode',[],'nBitCodes',[],'date_num',[],'timestamp',[],'switch_times',[],'switch_detected',[],'xyz_micron',[],'xyz_submicron',[],'piezo',[],'laser_power',[]);
        bitCodes=struct('nBitCodes',[],'scim_bitCodes_raw',[],'scim_bitCodes',[],'MWorks_bitCodes',[],'mwk_file_name','','event_codec',[],'offset',[],'max_val',[])
        FOV_info=struct('coords',[],'center',[],'Z_depth',[],'size_px',[],'pixel_size_micron',[],'size_um',[])
        
        
        motion_correction=struct('reference_image',struct('idx',[],'im',[],'shift_matrix',[],'min_val',[],'iBest',[],'total_shift',[]), ...
            'shift_matrix',[],'ignore_frames',[],'variability_threshold',[]);
        
        MIP_avg=struct('data',[],'gamma_val',[]);
        MIP_max=struct('data',[],'gamma_val',[]);
        MIP_std=struct('data',[],'gamma_val',[]);
        MIP_cc_local=struct('data',[],'gamma_val',[]);
        
        Experiment_info=struct('exp_type',[],'exp_name','','stimulus_data',struct('timestamp',[]));
        
        ROI_definitions=struct('ROI',struct(...
            'ROI_nr',[],'base_coord',[],'nCoords',[],'coords',[],...
            'ellipse_properties',[],'ellipse_coords',[],'coords_MIP',[],'coords_MIP_plot',[],'center_coords',[],'ellipse_coords_centered',[],...
            'ROI_rect',[],'mask_soma',[],'mask_neuropil',[],'timeseries_soma',[],'time_series_neuropil',[])...
            );
        
        Activity_traces=struct('activity_matrix',[],...
            'normalization_matrix',[],...
            'spike_matrix',[],...
            'oopsi_options',[],...
            'extraction_options',struct(...
            'neuropil_subtraction', struct('use',1,'factor',.7),...
            'drift_correction',     struct('use',1,'averaging_interval',15,'prctile',8),...
            'calc_delta_f_method',  6,...
            'do_FastNegDeconv',     1)...
            );
        
        Results=struct('RF_maps',[],'RSVP',[])
        
        Z_stack=struct('CC_matrix',[],'im',[],'est_Z_depth',[],'CC_val',[]);
        
        motion_proxy=[];
                
        green=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];
        
        last_action='';
        elapsed=[];
    end
    
    
    methods
        %%% Constructor
        function self=imaging_dataset(varargin)
            self.file_name=varargin{1};
            
            %%% split up in root_folder (until 2p-data)
            idx=strfind(self.file_name,self.folder_info.main_folder);
            if idx>0
                self.folder_info.root_folder=self.file_name(1:idx+7);
                self.folder_info.rel_path=self.file_name(idx+8:end);
            else % unable to parse out root folder
                self.folder_info.root_folder='';
                self.folder_info.rel_path=self.file_name;
            end
            [self.folder_info.data_folder,self.folder_info.raw_name]=fileparts(self.file_name);
            self.folder_info.save_folder=fullfile(self.folder_info.data_folder,self.folder_info.save_folder_root);
            self.save_name=fullfile(self.folder_info.save_folder,[self.folder_info.raw_name '.mat']);
        end
        
        %%% In case we need to change the root_folder to raw files
        function rebase(varargin)
            self=varargin{1};
            if nargin>=2
                self.folder_info.root_folder=varargin{2};
            else
                error('Rebase() needs folder to be used as root...')
            end
            self.file_name=fullfile(self.folder_info.root_folder,self.folder_info.rel_path);
        end
        
        
        %%% Extract info from movie
        function get_mov_info(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.mov_info.nFrames)
                %%% Get basic info about movie
                info=imfinfo(self.file_name); % never save info, is huge
                self.mov_info.nFrames=length(info);
                self.mov_info.Width=info(1).Width;
                self.mov_info.Height=info(1).Height;
            else
                disp('Using existing mov_info')
            end
            
            self.elapsed=toc;
            self.last_action='get_mov_info';
        end
        
        function get_scim_data(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.mov_info.state)
                info=imfinfo(self.file_name);
                
                %% Get SCIM headers
                scim_info=info(1).ImageDescription;
                scinfo=strsplit(scim_info,char(13));
                scim_info=strjoin(scinfo,[';' char(13)]);
                eval([scim_info  ';']); % revive scan image variables
                self.mov_info.state=state;
                self.mov_info.frame_rate=state.acq.frameRate;
                self.mov_info.mov_start_time=state.internal.softTriggerTimeString;
            else
                disp('Using existing scim_info...')
            end
            self.elapsed=toc;
            self.last_action='get_scim_data';
        end
        
        function read_flyback(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.frame_info(1).nBitCodes)
                info=imfinfo(self.file_name); % never save info, is huge
                
                N=[];
                for iFrame=1:self.mov_info.nFrames
                    % read out last line of the frame, here we store frame specific
                    % info: bitCodes, position, laser power
                    flyback_line=double(imread(self.file_name,iFrame,'info',info,'PixelRegion',{[self.mov_info.Height self.mov_info.Height],[1 self.mov_info.Width]}));
                    
                    [a,N]=parse_flyback_line(flyback_line,N);
                    self.frame_info(iFrame)=a;
                end
                self.bitCodes.nBitCodes=N;
            else
                disp('Using existing frame_info...')
            end
            self.elapsed=toc;
            self.last_action='read_flyback';
        end
        
        
        %%% Get FOV info
        function get_FOV_info(varargin)
            self=varargin{1};
            if nargin>=2
                self.FOV_info.pixel_size_micron=varargin{2};
            else
                self.FOV_info.pixel_size_micron=.85;
            end
            if self.is_static_FOV()
                self.FOV_info.coords=self.frame_info(1).xyz_submicron;
                self.FOV_info.center=self.FOV_info.coords(1:2);
                self.FOV_info.Z_depth=self.FOV_info.coords(3);
                self.FOV_info.size_px=[self.mov_info.Width self.mov_info.Height];
                self.FOV_info.size_um=self.FOV_info.size_px*self.FOV_info.pixel_size_micron;
            end
            %
        end
        
        function static_FOV=is_static_FOV(varargin)
            self=varargin{1};
            if nargin>=2
                TH=varargin{2};
            else
                TH=3;
            end
            M=cat(1,self.frame_info.xyz_micron);
            if any(std(M)>TH)
                static_FOV=0;
                fprintf('Some axis shows more variation (std XYZ=[%3.2f %3.2f %3.2f]) than excepted (TH=%3.2f)...\n',[std(M) TH])
            else % treat data as FOV recording where ROIs can be defined
                static_FOV=1;
            end
        end
        
        
        %%% BitCode stuff
        function join_bitCodes(varargin)
            self=varargin{1};
            %self.bitCodes.scim_bitCodes_raw=medfilt1(cat(1,self.frame_info.bitCode_vector),11);
            self.bitCodes.scim_bitCodes_raw=cat(1,self.frame_info.bitCode_vector); % med filt does not work if transition goes 1111244444, will still leave 2 in there, need something asymmetrical
            %all(eq(cat(1,self.frame_info.bitCode_vector),self.bitCodes.scim_bitCodes_raw))
        end
        
        function clean_up_bitCodes_raw(varargin)
            self=varargin{1};
            N=self.bitCodes.nBitCodes;
            B=mode(reshape(self.bitCodes.scim_bitCodes_raw,N,[]));
            out=repmat(B,N,1);
            self.bitCodes.scim_bitCodes_raw=out(:);
        end
                
        function get_scim_bitCodes(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.bitCodes.scim_bitCodes)
                self.join_bitCodes()
                self.clean_up_bitCodes_raw() % to get rid of transitions
                A=self.bitCodes.scim_bitCodes_raw;
                breaks=find([-1 ; diff(A)]);
                
                %                 sel=diff([breaks;0])<149;
                %
                %                 if length(sel)==length(breaks)
                %                     breaks(sel)=[];
                %                 else
                %                     error('messed up...')
                %                 end
                
                bitCode_time=(1/self.mov_info.frame_rate)/self.bitCodes.nBitCodes;
                T=((0:length(A)-1))'*bitCode_time;
                
                self.bitCodes.scim_bitCodes=[T(breaks) A(breaks) [0;diff(T(breaks))]];
            else
                disp('Using existing scim_bitCodes...')
            end
            
            self.elapsed=toc;
            self.last_action='get_scim_bitCodes';
        end
        
        function get_MWorks_bitCodes(varargin)
            tic
            self=varargin{1};
            
            if ismac
                if isempty(self.bitCodes.MWorks_bitCodes)
                    if nargin>=2
                        F=varargin{2};
                    else % if none specified, use first mwk file found in data_folder
                        %F='2015-07-17_AG02.mwk';
                        files=scandir(self.folder_info.data_folder,'.mwk');
                        F=files(1).name;
                    end
                    mwk_file_name=fullfile(self.folder_info.data_folder,F);
                    %disp('Reading MWK file...')
                    A=getCodecs(mwk_file_name);
                    event_codec=A.codec;
                    
                    %%% Get stim update events
                    tag_name='#stimDisplayUpdate';
                    [MW_events,nEvents]=get_events_by_name(mwk_file_name,tag_name,event_codec);
                    
                    MW_bitCodes=zeros(nEvents,3);
                    for iEvent=1:nEvents
                        event=MW_events(iEvent);
                        register=1;
                        if length(event.data)==3
                            bitCode=event.data{3}.bit_code;
                        elseif length(event.data)==2
                            bitCode=event.data{2}.bit_code;
                        else
                            % ignore
                            register=0;
                        end
                        if register==1
                            MW_bitCodes(iEvent,:)=[double(event.time_us)/1e6 double(bitCode) [0;diff(double(event.time_us)/1e6)]];
                        end
                    end
                    
                    self.bitCodes.MWorks_bitCodes=MW_bitCodes;
                    self.bitCodes.mwk_file_name=mwk_file_name;
                    self.bitCodes.event_codec=event_codec;
                else
                    disp('Using existing MWorks_bitCodes...')
                end
            else
                disp('MWorks magic only works in Matlab for mac...')
            end
            self.elapsed=toc;
            self.last_action='get_MWorks_bitCodes';
        end
        
        function find_offset(varargin)
            tic
            self=varargin{1};
            if isempty(self.bitCodes.offset)
                A=self.bitCodes.scim_bitCodes(:,2);
                B=self.bitCodes.MWorks_bitCodes(:,2);
                CC=normxcorr2(A,B);
                [self.bitCodes.max_val,loc]=max(CC);
                if self.bitCodes.max_val>.99
                    self.bitCodes.offset=loc-length(A)+1;
                else
                    self.bitCodes.offset=[];
                end
            else
                disp('Using existing offset...')
            end
            
            self.elapsed=toc;
            self.last_action='find_offset';
        end
        
        
        %%% Methods related to the raw movie
        function frames=get_frames(varargin)
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                idx=varargin{2};
            else
                idx=1:self.mov_info.nFrames;
            end
            if nargin>=3&&~isempty(varargin{3})
                apply_motion_correction=varargin{3};
            else
                apply_motion_correction=0;
            end
            N=length(idx);
            info=imfinfo(self.file_name); % never save info, is huge
            frames=zeros([self.mov_info.Height-1 self.mov_info.Width N]);
            for iFrame=1:N
                frame=double(imread(self.file_name,idx(iFrame),'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                if apply_motion_correction==1
                    offset_ij=-self.motion_correction.shift_matrix(iFrame,[2 3]);
                    frame=offsetIm(frame,offset_ij(1),offset_ij(2),0);
                else
                    % do nothing
                end
                if size(frame,3)==3
                    frame=rgb2gray(frame/256);
                    %frame=frame(:,:,1);
                end
                frames(:,:,iFrame)=frame;
            end
        end
       
        function find_blank_frames(varargin)
            %%% Check for unusually dark frames, laser power not turned up
            %%% yet.
            tic
            self=varargin{1};
            if isempty(self.mov_info.mean_lum)
                info=imfinfo(self.file_name); % never save info, is huge
                
                self.mov_info.mean_lum=zeros(self.mov_info.nFrames,1);
                for iFrame=1:self.mov_info.nFrames
                    frame=double(imread(self.file_name,iFrame,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                    self.mov_info.mean_lum(iFrame,1)=mean(frame(:));
                end
                
                self.mov_info.blank_frames=false(size(self.mov_info.mean_lum));
                
                % Search for longer periods that might cause previous condition to fail
                % low std, low mean
                z_scores=(self.mov_info.mean_lum(1:round(end/4))-mean(self.mov_info.mean_lum(round(end/4):end)))/std(self.mov_info.mean_lum(round(end/4):end));
                self.mov_info.blank_frames(1:find(z_scores<-4,1,'last')+1)=true;
            else
                disp('Using existing blank_frames...')
            end
            
            self.elapsed=toc;
            self.last_action='find_blank_frames';
        end
        
        
        %%% Get experiment info from MWorks
        function get_exp_type(varargin)
            self=varargin{1};
            
            
            if ismac
                exp_type=[];
                exp_name='';
                
                %%% Try first to read from the MWorks variable exptype directly
                tag_name='ExpType';
                expType_events=get_events_by_name(self.bitCodes.mwk_file_name,tag_name,self.bitCodes.event_codec);
                if ~isempty(expType_events)
                    exp_name_vector={'RSVP','Retinomapping'};
                    exp_type=mode(expType_events);
                    exp_name=exp_name_vector{exp_type};
                end
                
                %%% Fallback is to look for the existance of specific tags
                %%% unique to either experiment
                tag_names={self.bitCodes.event_codec.tagname}';
                tag_nr=find(ismember(tag_names,'stm_pos_x'),1);
                if ~isempty(tag_nr)
                    exp_type=1; % RSVP
                    exp_name='RSVP';
                end
                tag_nr=find(ismember(tag_names,'show_vertical_bar'),1);
                if ~isempty(tag_nr)
                    exp_type=2; % Retinomapping
                    exp_name='Retinomapping';
                end
                
                if isempty(exp_type)
                    disp('Unable to determine experiment type...')
                end
                
                self.Experiment_info.exp_type=exp_type;
                self.Experiment_info.exp_name=exp_name;
            else
                % no option
            end
        end
        
        function get_MWorks_stimulus_info(varargin)
            self=varargin{1};
            
            if mean(eq(self.bitCodes.scim_bitCodes(:,2),self.bitCodes.MWorks_bitCodes(self.bitCodes.offset:self.bitCodes.offset-1+size(self.bitCodes.scim_bitCodes,1),2)))>.99
                stim_times=self.bitCodes.MWorks_bitCodes(self.bitCodes.offset:self.bitCodes.offset-1+size(self.bitCodes.scim_bitCodes,1),1);
                stim_times([1 end]); % these time are sufficient to capture all events for this experiment
                
                if ~isempty(self.Experiment_info.exp_type)
                    switch self.Experiment_info.exp_type
                        case 1
                            %%% Read all events between start and end of
                            %%% session
                            tag_name='#stimDisplayUpdate';
                            [MW_events,nEvents]=get_events_by_name(self.bitCodes.mwk_file_name,tag_name,self.bitCodes.event_codec,uint64(stim_times(1)*1e6),uint64(stim_times(end)*1e6));
                            
                            %stimulus_data=self.Experiment_info.stimulus_data;
                            for iEvent=1:nEvents
                                event=MW_events(iEvent);
                                stim_duration=diff([self.bitCodes.MWorks_bitCodes(self.bitCodes.offset-1+iEvent) self.bitCodes.MWorks_bitCodes(self.bitCodes.offset+iEvent)]);
                                
                                stimulus_data(iEvent).timestamp=event.time_us;
                                stimulus_data(iEvent).bit_code=event.data{end}.bit_code;
                                stimulus_data(iEvent).stim_duration=stim_duration;
                                switch length(event.data)
                                    case 2 % blank
                                        stimulus_data(iEvent).stim_present=0;
                                    case 3 % stim
                                        stimulus_data(iEvent).stim_present=1;
                                        stimulus_data(iEvent).stim_id=str2double(event.data{2}.name)+1;
                                        stimulus_data(iEvent).size_x=event.data{2}.size_x;
                                        stimulus_data(iEvent).size_y=event.data{2}.size_y;
                                        stimulus_data(iEvent).pos_x=event.data{2}.pos_x;
                                        stimulus_data(iEvent).pos_y=event.data{2}.pos_y;
                                        stimulus_data(iEvent).position_nr=self.get_position_number([stimulus_data(iEvent).pos_x stimulus_data(iEvent).pos_y]);
                                        stimulus_data(iEvent).rotation=event.data{2}.rotation;
                                        nPositions=32;
                                        stimulus_data(iEvent).condition_nr=(stimulus_data(iEvent).stim_id-1)*nPositions+stimulus_data(iEvent).position_nr;
                                    otherwise
                                        error('Unexpected amount of fields in event data...')
                                end
                            end
                            self.Experiment_info.stimulus_data=stimulus_data;
                            
                        case 2 % Retinomapping
                            error('No processing pipeline defined for this experiment')
                        otherwise
                            error('No processing pipeline defined for this experiment')
                    end
                else
                    error('Experiment type not defined...')
                end
                
            else
                error('Bit code mismatch...')
            end
        end
        
        function position_nr=get_position_number(varargin)
            %self=varargin{1};
            coords=varargin{2};
            positions_LUT=[-45.5,-19.5;-45.5,-6.5;-45.5,6.5;-45.5,19.5;-32.5,-19.5;-32.5,-6.5;-32.5,6.5;-32.5,19.5;-19.5,-19.5;-19.5,-6.5;-19.5,6.5;-19.5,19.5;-6.5,-19.5;-6.5,-6.5;-6.5,6.5;-6.5,19.5;6.5,-19.5;6.5,-6.5;6.5,6.5;6.5,19.5;19.5,-19.5;19.5,-6.5;19.5,6.5;19.5,19.5;32.5,-19.5;32.5,-6.5;32.5,6.5;32.5,19.5;45.5,-19.5;45.5,-6.5;45.5,6.5;45.5,19.5];
            position_nr=find(ismember(positions_LUT,coords,'rows'));
        end
        
        function create_stim_matrix(varargin)
            %%% This method will take us from a structure with events to a
            %%% matrix with stimulus properties one line for each frame
            %%% recorded.
            self=varargin{1};
            stimulus_data=self.Experiment_info.stimulus_data;
            % map each event to a range of frames
            %A=self.bitCodes.scim_bitCodes; % fixed!
            B=mode(reshape(self.bitCodes.scim_bitCodes_raw,50,[]))';
            
            trial_mapping=self.expand_trial_numbers(B);
            stim_matrix=zeros(self.mov_info.nFrames,10)-1;
            for iTrial=1:length(stimulus_data)
                sel=trial_mapping==iTrial;
                stim_props=[iTrial stimulus_data(iTrial).stim_present stimulus_data(iTrial).condition_nr stimulus_data(iTrial).stim_id stimulus_data(iTrial).position_nr stimulus_data(iTrial).pos_x stimulus_data(iTrial).pos_y stimulus_data(iTrial).size_x stimulus_data(iTrial).size_y stimulus_data(iTrial).rotation];
                stim_matrix(sel,1:length(stim_props))=repmat(stim_props,sum(sel),1);
            end
            stim_matrix=[(1:self.mov_info.nFrames)' stim_matrix];
            self.Experiment_info.stim_matrix=stim_matrix;
        end
        
        function res=expand_trial_numbers(varargin)
            self=varargin{1};
            vector=varargin{2};
            cur_val=[];
            count=0;
            res=zeros(size(vector));
            for index=1:length(vector)
                val=vector(index);
                if val==cur_val
                    % use same trial number
                else % increase trial number
                    cur_val=val;
                    count=count+1;
                end
                res(index)=count;
            end
        end
                
        
        %%%% Motion correction methods
        function do_motion_detection(varargin)
            %%% Implement frame to frame difference as a proxy for animal
            %%% motion
            % Not used often
            tic
            self=varargin{1};
            info=imfinfo(self.file_name); % never save info, is huge
            
            diff_vector=zeros(self.mov_info.nFrames,1);
            fprintf('Completed: %03d%%\n',0)
            for iFrame=1:self.mov_info.nFrames-1
                cur_frame=double(imread(self.file_name,iFrame,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                next_frame=double(imread(self.file_name,iFrame+1,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                
                diff_vector(iFrame+1)=sum(sqrt((cur_frame(:)-next_frame(:)).^2));
                del_str=repmat('\b',1,5);
                fprintf([del_str '%03d%%\n'],round(iFrame/(self.mov_info.nFrames-1)*100))
            end
            self.motion_proxy=diff_vector;
            
            self.elapsed=toc;
            self.last_action='do_motion_correction';
        end
        
        function find_reference_image(varargin)
            % Implementing motion correction method by Poort et al. 2015 Neuron
            tic
            self=varargin{1};
            
            nBlocks=30;
            block_size=10; % less frames, but more time given our lower sampling rate: 3 vs 30Hz
            nSamples=100;
            
            if isempty(self.motion_correction.reference_image.im)
                %%% Get frames
                frames=self.get_frames();
                
                start_frame=find(self.mov_info.blank_frames==0,1,'first'); % find first non-blank frame
                end_frame=self.mov_info.nFrames;
                
                %%% space 30 blocks of 30 linearly across the stack (or 1s)
                A=floor(linspace(start_frame,end_frame-block_size,nBlocks));
                M=[A(:) A(:)+block_size];
                
                fprintf('Calculating reference image... ')
                reference_image_candidates=zeros([size(frames,1) size(frames,2) nBlocks]);
                for iBlock=1:nBlocks
                    idx=M(iBlock,1):M(iBlock,2);
                    reference_image_candidates(:,:,iBlock)=mean(frames(:,:,idx),3);
                end
                
                % space 100 single frames linearly across the stack
                A=floor(linspace(start_frame,end_frame,nSamples));
                sample_images=frames(:,:,A);
                
                % cross-correlate
                %self.motion_correction.reference_image.CC_matrix=zeros(nSamples,nBlocks);
                shift_matrix=zeros(nSamples,2,nBlocks);
                for iBlock=1:nBlocks
                    ref=reference_image_candidates(:,:,iBlock);
                    for iSample=1:nSamples
                        sample=sample_images(:,:,iSample);
                        % smooth B?
                        [r,c]=PCdemo(ref,sample);
                        shift_matrix(iSample,:,iBlock)=[c r];
                    end
                end
                
                % find best
                total_shift=zeros(nBlocks,1);
                for iBlock=1:nBlocks
                    total_shift(iBlock)=sum(sqrt(sum(shift_matrix(:,:,iBlock).^2,2)));
                end
                [min_val,iBest]=min(total_shift);
                
                self.motion_correction.reference_image.idx=M(iBest,1):M(iBest,2);
                self.motion_correction.reference_image.im=reference_image_candidates(:,:,iBest);
                self.motion_correction.reference_image.shift_matrix=shift_matrix;
                self.motion_correction.reference_image.min_val=min_val;
                self.motion_correction.reference_image.iBest=iBest;
                self.motion_correction.reference_image.total_shift=total_shift;
                
                self.elapsed=toc;
                self.last_action='find_reference_image';
                
                fprintf('Done!\n')
            else
                disp('Using existing reference image...')
            end
        end
        
        function do_motion_correction(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.motion_correction.shift_matrix)
                info=imfinfo(self.file_name); % never save info, is huge
                
                ref=self.motion_correction.reference_image.im;
                
                fprintf('Running motion correction... ')
                str=sprintf('%03d%%\n',0);
                fprintf(str)
                shift_matrix=zeros(self.mov_info.nFrames,3);
                for iFrame=1:self.mov_info.nFrames
                    cur_frame=double(imread(self.file_name,iFrame,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                    if size(cur_frame,3)==3
                        cur_frame=rgb2gray(cur_frame/256);
                    end
                    
                    [r,c]=PCdemo(ref,cur_frame);
                    shift_matrix(iFrame,:)=[iFrame c r];
                    
                    del_str=repmat('\b',1,length(str));
                    fprintf([del_str '%03d%%\n'],round(iFrame/(self.mov_info.nFrames-1)*100))
                end
                self.motion_correction.shift_matrix=shift_matrix;
            else
                disp('Using existing shift_matrix...')
            end
            
            self.elapsed=toc;
            self.last_action='do_motion_correction';
        end
        
        function reset_motion_correction(varargin)
            self=varargin{1};
            self.motion_correction.shift_matrix=[];
        end
        
        function find_motion_frames(varargin)
            tic
            self=varargin{1};
            if nargin>=2
                TH=varargin{2};
            else
                TH=10;
            end
            
            %T=self.motion_correction.shift_matrix(:,1);
            D=calc_dist([self.motion_correction.shift_matrix(:,2:3)*0 self.motion_correction.shift_matrix(:,2:3)]);
            M=[cat(1,D,[0; 0;0]) cat(1,0,D,[0;0]),cat(1,[0; 0],D,0),cat(1,[0; 0;0],D)];M(end-2:end,:)=[];
            self.motion_correction.ignore_frames=std(M,[],2)>TH;
            self.motion_correction.variability_threshold=TH;
            
            self.elapsed=toc;
            self.last_action='find_motion_frames';
        end
        
        
        %%% Calc MIPs using motion correction
        function do_calc_MIPs(varargin)
            tic
            self=varargin{1};
            
            debugging=1;
            
            %if any([isempty(self.MIP_avg.data) isempty(self.MIP_max.data) isempty(self.MIP_std.data) isempty(self.MIP_cc_local.data)])
            if any([isempty(self.MIP_avg.data) isempty(self.MIP_max.data) isempty(self.MIP_std.data)])
                fprintf('Calculating MIPs... ')
                %%% Get motion corrected frames
                frames=self.get_frames([],1);
                
                %%% Remove bad frames
                del=self.mov_info.blank_frames==1|self.motion_correction.ignore_frames==1;
                frames(:,:,del)=[];
                
                %%% Do median filter to get rid of speckles
                if size(frames,3)<2000
                    frames=medfilt3(frames,[1 1 3]);
                end
                
                %%% Calc all MIPs
                if isempty(self.MIP_avg.data)
                    self.MIP_avg.data=mean(frames,3);
                    self.MIP_avg.gamma_val=1;
                end
                if isempty(self.MIP_max.data)
                    self.MIP_max.data=max(frames,[],3);
                    self.MIP_max.gamma_val=.6;
                end
                if isempty(self.MIP_std.data)
                    self.MIP_std.data=std(frames,[],3);
                    self.MIP_std.gamma_val=.8;
                end
                if isempty(self.MIP_cc_local.data)
                    if debugging==0 % only if we are not debugging, takes a long time to run
                        self.MIP_cc_local.data=CrossCorrImage(frames);
                        self.MIP_cc_local.gamma_val=.6;
                    end
                end
                fprintf('Done!\n')
            else
                disp('Using existing MIPs...')
            end
            
            self.elapsed=toc;
            self.last_action='do_calc_MIPs';
        end
        
        
        %%% Do ROI extraction
        function do_trace_extraction(varargin)
            tic
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                ROI_definition_nr=varargin{2};
            else
                ROI_definition_nr=1;
            end
            
            if isempty(self.Activity_traces.activity_matrix)
                %% Grab average ROI activity over the motion corrected movie
                %if ~isfield(self.ROI_definitions(ROI_definition_nr).ROI,'timeseries_soma')
                temp=get_ROI_definitions(self,ROI_definition_nr);
                if isempty(temp)
                    disp('No valid ROI definitions found...')
                else
                    if ~isfield(temp(end),'timeseries_soma')||isempty(temp(end).timeseries_soma)
                        ROIs=self.extract_traces(ROI_definition_nr);
                        self.ROI_definitions(ROI_definition_nr).ROI=ROIs;
                    else % load existing definitions
                        ROIs=get_ROI_definitions(self,ROI_definition_nr);
                    end
                    %ROI_vector=cat(1,ROIs.ROI_nr);
                    nROI=length(ROIs);
                    nFrames=self.mov_info.nFrames;
                    
                    activity_matrix=zeros(nROI,nFrames);
                    spike_matrix=activity_matrix;
                    for iROI=1:nROI
                        trace=self.process_trace(ROIs,iROI);
                        %subplot(nROI,1,iROI)
                        %plot(trace)
                        
                        activity_matrix(iROI,:)=trace;
                        
                        if self.Activity_traces.extraction_options.do_FastNegDeconv==1
                            %if self.Activity_traces.do_FastNegDeconv==1
                            %% Fast oopsi
                            oopsi_options.Ncells=1;
                            oopsi_options.T=nFrames;
                            oopsi_options.dt=1/self.mov_info.frame_rate*0.7;
                            [n_best, P_best, oopsi_options, C]=fast_oopsi(trace,oopsi_options);
                            spike_matrix(iROI,:)=n_best;
                        end
                    end
                    
                    self.Activity_traces.activity_matrix=activity_matrix';
                    self.Activity_traces.spike_matrix=spike_matrix';
                    self.Activity_traces.oopsi_options=oopsi_options;
                    
                    self.elapsed=toc;
                    self.last_action='do_trace_extraction';                    
                end
            else
                disp('Using existing trace_matrix')
            end
        end
        
        function reset_trace_matrix(varargin)
            self=varargin{1};
            self.Activity_traces.activity_matrix=[];
        end
            
        function ROIs=extract_traces(varargin)
            %%% Goes way faster than before!
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                ROI_definition_nr=varargin{2};
            else
                ROI_definition_nr=1;
            end
            
            tic
            fprintf('Loading frames... ')
            frames=self.get_frames([],1); % get all frames, apply motion correction
            %nFrames=size(frames,3);
            fprintf('took: %3.1fs\n',toc);
            
            ROIs=get_ROI_definitions(self,ROI_definition_nr);
            nROI=length(ROIs);
            tic
            fprintf('Extracting %d ROIs data... ',nROI)
            for iROI=1:nROI
                rect=ROIs(iROI).ROI_rect;
                vol=frames(rect(2):rect(4),rect(1):rect(3),:);
                mask=repmat(ROIs(iROI).mask_soma,1,1,size(vol,3));
                res=vol.*mask;
                ROIs(iROI).timeseries_soma=squeeze(mean(mean(res,1),2));
                
                mask=repmat(ROIs(iROI).mask_neuropil,1,1,size(vol,3));
                res=vol.*mask;
                ROIs(iROI).timeseries_neuropil=squeeze(mean(mean(res,1),2));
            end
            fprintf('took: %3.1fs\n',toc);
        end
        
        function trace=process_trace(varargin)
            self=varargin{1};
            ROIs=varargin{2};
            iROI=varargin{3};
            ROI=ROIs(iROI);
            
            extraction_options=self.Activity_traces.extraction_options;
            frame_rate=self.mov_info.frame_rate;
            nFrames=self.mov_info.nFrames;
            blank_frames=self.mov_info.blank_frames;
            if extraction_options.neuropil_subtraction.use
                %%% Step 01: get F
                F_raw=ROI.timeseries_soma;
                F_neuropil=ROI.timeseries_neuropil;
                
                if any(blank_frames)
                    F_raw(blank_frames)=mean(F_raw(~blank_frames));
                    F_neuropil(blank_frames)=mean(F_neuropil(~blank_frames));
                end
                
                %%% subtract neuropil signals
                F=F_raw-F_neuropil*extraction_options.neuropil_subtraction.factor; % values .5 (Feinberg) and .7 (Kerlin) have been reported: .7 is used for 16x Nikon
            else
                %error('No neuropil subtraction?')
                F=ROI.timeseries_soma;
            end
            
            if extraction_options.drift_correction.use
                nFrames_avg=extraction_options.drift_correction.averaging_interval*round(frame_rate);
                
                % make sure this number even
                if rem(nFrames_avg,2)==1
                    nFrames_avg=nFrames_avg+1;
                end
                
                % take nFrames_avg around each sample and get 8th prctile value
                % of these values.
                prctile_vector=zeros(size(F));
                for iSample=1:nFrames
                    if between(iSample,[nFrames_avg/2+1 nFrames-nFrames_avg/2])
                        hist_values=F(iSample-nFrames_avg/2+1:iSample+nFrames_avg/2);
                    elseif iSample<nFrames_avg % handles start of trace
                        hist_values=F(1:iSample+nFrames_avg/2);
                    elseif iSample>nFrames-nFrames_avg % handles end of trace
                        hist_values=F(iSample-nFrames_avg/2+1:end);
                    end
                    prctile_vector(iSample)=prctile(hist_values,extraction_options.drift_correction.prctile);
                end
                % subtract prctile vector from signal
                F_no_drift=F-prctile_vector;
            else
                error('No drift correction?')
            end
            
            
            %switch self.Activity_traces.calc_delta_f_method
            switch extraction_options.calc_delta_f_method
                case 6
                    % use rolling average to select no-modulation parts
                    F_no_drift_smooth=smooth(F_no_drift(blank_frames==0),round(frame_rate*4));
                    
                    TH=prctile(F_no_drift_smooth,2)+512; % optimize!!!
                    selection_vector=F_no_drift_smooth<TH;
                    mu=mean(F_no_drift(selection_vector));
                    sigma=std(F_no_drift(selection_vector));
                    
                    % Store normalizaion values
                    self.Activity_traces.normalization_matrix(iROI,1:4)=[iROI extraction_options.calc_delta_f_method mu sigma];
                    
                    delta_F_no_drift=(F_no_drift-mu)/mu;
                    y_label='\DeltaF/F';
                    fixed_y_scale=30;
                otherwise
                    error('Still need to port other methods...')
            end
            extraction_options.y_label=y_label;
            extraction_options.fixed_y_scale=fixed_y_scale;
            %self.Activity_traces.normalization_matrix=normalization_matrix;
            
            delta_F_no_drift(blank_frames)=0;
            
            trace=delta_F_no_drift;            
        end                
        
        
        %%% Combine activity and stim properties
        function combine_act_stim(varargin)
            self=varargin{1};
            if nargin>=2
                iROI=varargin{2};
            else
                iROI=1;
            end
            if nargin>=3
                variable=varargin{3};
            else
                variable=5; % by stim
            end
            
            Stim=self.Experiment_info.stim_matrix;
            Resp=self.Activity_traces.activity_matrix;
            %Resp=self.Activity_traces.spike_matrix;
                        
            %plot(self.mov_info.mean_lum)
            %plot([zscore(Stim(:,5)) zscore(self.mov_info.mean_lum)])
            %corr([zscore(Stim(:,5)) zscore(self.mov_info.mean_lum)])
            
            
            %self.analyse_data([Stim(:,1:3) Stim(:,variable) self.mov_info.mean_lum])
            %self.analyse_data_new([Stim(:,1:3) Stim(:,variable) self.mov_info.mean_lum])
            
            %self.analyse_data([Stim(:,1:3) Stim(:,variable) Resp(:,iROI)])
            self.analyse_data_new([Stim(:,1:3) Stim(:,variable) Resp(:,iROI)])            
        end
        
        function analyse_data_new(varargin)
            self=varargin{1};
            data=varargin{2};
            condition_vector=data(:,4);
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
                condition_matrix(iTrial,end)=mean(data(frame_range,end));
                
                frame_range=trial_info(2)+frame_selector_trace(1):trial_info(2)+frame_selector_trace(2);
                calcium_matrix(iTrial,:)=data(frame_range,end);
            end
            M=pivotTable(condition_matrix,5,'mean',7);
            S=pivotTable(condition_matrix,5,'std',7);
            E=pivotTable(condition_matrix,5,'ste',7);
            %[M S E]
            
            figure(2)
            subplot(211)
            plot(data(:,end))
            subplot(212)
            if nConditions==32
                MAP=flipud(reshape(M,4,8));
                imagesc(MAP)
                set(gca,'CLim',[-3 3])
                self.Results.RF_maps(:,:,iROI)=MAP;
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
            end
            
        end
        
        function analyse_data(varargin)
            self=varargin{1};
            M=varargin{2};
            
            %%% Average data per presentation
            condition_vector=M(M(:,3)==1,4);
            %parse_conditions(condition_vector) % how did we use this?
            conditions=unique(condition_vector);
            nConditions=length(conditions);
            condition_matrix=zeros(nConditions,5);
            condition_matrix_detail=struct;
            for iCondition=1:nConditions
                condition_nr=conditions(iCondition);
                
                % find all repeats for this condition
                trial_numbers=unique(M(M(:,4)==condition_nr,2));
                nRepeats=length(trial_numbers);
                
                repeat_matrix=zeros(nRepeats,3);
                for iRepeat=1:nRepeats
                    trial_nr=trial_numbers(iRepeat);
                    sel=M(:,2)==trial_nr;
                    frame_numbers=M(sel,1);
                    
                    switch 1
                        case 1
                            resp_data=M(sel,end);
                        case 2
                            start_frame=frame_numbers(1);
                            end_frame=frame_numbers(end);
                            frame_idx=start_frame+2:end_frame+4;
                            resp_data=M(frame_idx,end);
                    end
                    repeat_matrix(iRepeat,:)=[iRepeat mean(resp_data) std(resp_data)];
                end
                condition_matrix(iCondition,:)=[iCondition nRepeats mean(repeat_matrix(:,2)) std(repeat_matrix(:,2)) ste(repeat_matrix(:,2))];
                condition_matrix_detail(iCondition).repeat_matrix=repeat_matrix;
            end
            
            %%% Bar plot
            figure(2)
            subplot(211)
            plot(M(:,end))
            subplot(212)
            bar(condition_matrix(:,3))
            hold on
            errorbar(condition_matrix(:,3),condition_matrix(:,5),'r.')
            hold off
            if nConditions==32
                imagesc(flipud(reshape(condition_matrix(:,3),4,8)))
            end
            
            %%% ANOVA
            %[P,T,STATS,TERMS]=anovan(M(:,2),M(:,1))
        end
        
        
        %%% File i/o
        function save_data(varargin)
            self=varargin{1};
            
            if ~isdir(self.folder_info.save_folder)
                mkdir(self.folder_info.save_folder)
            end
            session_data=self;
            save(self.save_name,'session_data')
            
            fprintf('Data saved to %s!\n',self.save_name)
        end
        
        function load_data(varargin)
            self=varargin{1};
            
            load(self.save_name,'session_data')
            self=session_data;
            
            fprintf('Data loaded from %s!\n',self.save_name)
        end
        
        
        %%% Visualization
        function imshow(varargin)
            self=varargin{1};
            im=varargin{2};
            if isstruct(im)
                gamma_val=im.gamma_val;
                im=im.data;
            else
                if nargin>=3
                    gamma_val=varargin{3};
                else
                    gamma_val=.5;
                end
            end
            imshow(calc_gamma(im,gamma_val),[])
            colormap(self.green)
        end
        
        function visualize_motion_correction(varargin)
            self=varargin{1};
            if nargin>=2
                apply_motion_correction=varargin{2};
            else
                apply_motion_correction=0;
            end
            frames=self.get_frames([],apply_motion_correction);
            ref=self.motion_correction.reference_image.im;
            center=[self.mov_info.Width self.mov_info.Height]/2;
            
            figure()
            H=imshow(calc_gamma(ref,.5),[]);
            hold on
            p=plot(center(1),center(2),'m*');
            
            colormap(self.green)
            hold off
            for iFrame=1:self.mov_info.nFrames
                cur_frame=frames(:,:,iFrame);
                coords=center+self.motion_correction.shift_matrix(iFrame,2:3);
                
                set(H,'Cdata',calc_gamma(cur_frame,.5))
                set(p,'xData',coords(1),'yData',coords(2))
                if self.motion_correction.ignore_frames(iFrame)==1
                    set(p,'marker','o','color','y')
                else
                    set(p,'marker','*','color','m')
                end
                drawnow
            end
        end
        
        function plot_traces(varargin)
            self=varargin{1};
            A=self.Activity_traces.activity_matrix;
            nROI=size(A,2);
            nCols=ceil(sqrt(nROI));
            nRows=ceil(nROI/nCols);
            figure()
            for iROI=1:nROI
                subplot(nRows,nCols,iROI)
                plot(A(:,iROI))
                axis([1 size(A,1) -10 40])
                set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            end
        end
        
        
        %%% Extra features
        function find_FOV_in_stack(varargin)
            tic
            self=varargin{1};
            stack=varargin{2};
            
            FOV=self.MIP_avg.data;
            
            if isempty(self.Z_stack.CC_matrix)
                frames=stack.get_frames(0); % no motion correction
                coords=cat(1,stack.frame_info.xyz_submicron);
                Z_depth=coords(:,3);
                
                dataMatrix=zeros(stack.mov_info.nFrames,6);
                t0=clock;
                for iFrame=1:stack.mov_info.nFrames
                    frame=frames(:,:,iFrame);
                    [CC_max,offset]=im_align(frame,FOV);
                    dataMatrix(iFrame,:)=[iFrame Z_depth(iFrame) mean(frame(:)) CC_max offset];
                    progress(iFrame,stack.mov_info.nFrames,t0)
                end
                plot(dataMatrix(:,3))
                
                self.Z_stack.CC_matrix=dataMatrix;
            else
                disp('Using existing CC_matrix...')
            end
            
            M=self.Z_stack.CC_matrix;
            [max_val,loc]=max(M(:,4));
            frames=stack.get_frames(loc-4:loc+5);
            frames=medfilt3(frames,[1 1 3]);
            self.Z_stack.im=mean(frames,3);
            self.Z_stack.est_Z_depth=M(loc,2);
            self.Z_stack.CC_val=max_val;
            
            self.elapsed=toc;
            self.last_action='find_FOV_in_stack';
            
        end
        
        function export_movie(varargin)
            self=varargin{1};
        end        
        
        function loop_status(varargin)
            self=varargin{1};
            
        end
    end
    
end