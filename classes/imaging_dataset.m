classdef imaging_dataset < handle
    properties
        file_name='';
        save_name='';
        folder_info=struct('main_folder','2p-data','root_folder','','rel_path','','data_folder','','raw_name','','save_folder','','save_folder_root','data_analysis')
        mov_info=struct('nFrames',[],'Width',[],'Height',[],'state',[],'frame_rate',[],'mov_start_time',[],'mean_lum',[],'blank_frames',[])
        frame_info=struct('version',[],'nBitCodes',[],'bitCode_vector',[],'main_bitCode',[],'date_num',[],'timestamp',[],'switch_times',[],'switch_detected',[],'xyz_micron',[],'xyz_submicron',[],'piezo',[],'laser_power',[]);
        bitCodes=struct('nBitCodes',[],'scim_bitCodes_raw',[],'scim_bitCodes',[],'MWorks_bitCodes',[],'mwk_file_name','','event_codec',[],'offset',[],'max_val',[])
        FOV_info=struct('coords',[],'center',[],'Z_depth',[],'size_px',[],'pixel_size_micron',[],'pixel_aspect_ratio',[],'size_um',[])
        
        
        motion_correction=struct('reference_image',struct('idx',[],'im',[],'shift_matrix',[],'min_val',[],'iBest',[],'total_shift',[]), ...
            'kernel',[],'shift_matrix',[],'motion_frames',[],'mismatch_frames',[],'ignore_frames',[],'variability_threshold',[]);
        
        MIP_avg=struct('data',[],'gamma_val',[]);
        MIP_max=struct('data',[],'gamma_val',[]);
        MIP_std=struct('data',[],'gamma_val',[]);
        MIP_cc_local=struct('data',[],'gamma_val',[]);
        MIP_kurtosis=struct('data',[],'gamma_val',[]);
        
        Experiment_info=struct('exp_type',[],'exp_name','','stim_duration',[],'stimulus_data',struct('timestamp',[]),'stim_matrix',[]);
        
        ROI_definitions=struct('ROI',struct(...
            'ROI_nr',[],'base_coord',[],'nCoords',[],'coords',[],...
            'ellipse_properties',[],'ellipse_coords',[],'coords_MIP',[],'coords_MIP_plot',[],'center_coords',[],'ellipse_coords_centered',[],...
            'ROI_rect',[],'mask_soma',[],'mask_neuropil',[],'timeseries_soma',[])...
            );
        ROI_definition_nr=[];
        
        ROI_matrix=struct;
        
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
        updated=0;
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
            self.updated=1;
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
            [self.folder_info.data_folder,self.folder_info.raw_name]=fileparts(self.file_name);
            self.folder_info.save_folder=fullfile(self.folder_info.data_folder,self.folder_info.save_folder_root);
            self.save_name=fullfile(self.folder_info.save_folder,[self.folder_info.raw_name '.mat']);
        end
        
        %%% In case we need to change the root_folder to raw files
        function rebase_tif(varargin)
            self=varargin{1};
            if nargin>=2
                new_folder=varargin{2};
            else
                error('Rebase_tif() needs folder to be used as root...')
            end
            self.file_name=fullfile(new_folder,[self.folder_info.raw_name '.tif']);
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
                
                self.elapsed=toc;
                self.last_action='get_mov_info';
                self.updated=1;
            else
                disp('Using existing mov_info')
            end
        end
        
        function get_scim_data(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.mov_info.state)
                info=imfinfo(self.file_name);
                
                %% Get SCIM headers
                scim_info=info(1).ImageDescription;
                if ~ismember(double(scim_info(end)),[10 13])
                    scim_info=[scim_info char(13)];
                    disp('Corrected last char')
                end
                scim_info=strrep(scim_info,char(10),char(13));
                scim_info=strrep(scim_info,char(13),[';' char(13)]);
                eval(scim_info); % revive scan image variables
                self.mov_info.state=state;
                
                self.mov_info.frame_rate=state.acq.frameRate;
                self.mov_info.mov_start_time=state.internal.softTriggerTimeString;
                
                self.elapsed=toc;
                self.last_action='get_scim_data';
                self.updated=1;
            else
                disp('Using existing scim_info...')
            end
            
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
                    try
                        self.frame_info(iFrame)=a;
                    catch
                        self.frame_info(iFrame)
                        a
                        die
                    end
                end
                
                N=mode(cat(1,self.frame_info.nBitCodes));
                
                for iFrame=1:self.mov_info.nFrames
                    V=self.frame_info(iFrame).bitCode_vector;
                    if length(V)<N
                        self.frame_info(iFrame).bitCode_vector=ones(N,1)*mode(V);
                        self.frame_info(iFrame).nBitCodes=length(self.frame_info(iFrame).bitCode_vector);
                        fprintf('Padded bitCodes in frame %d...\n', iFrame)
                    elseif length(V)>N
                        self.frame_info(iFrame).bitCode_vector=V(1:N);
                        self.frame_info(iFrame).nBitCodes=N;
                        fprintf('Cropped bitCodes in frame %d...\n', iFrame)
                    end
                end
                
                self.bitCodes.nBitCodes=N;
                
                self.elapsed=toc;
                self.last_action='read_flyback';
                self.updated=1;
            else
                disp('Using existing frame_info...')
            end
            
        end
        
        function recover_scim_bitCodes_from_MWorks(varargin)
            % for the rare cases flyback line gets overwritten, reconstruct
            % frame bitcodes from MWork timestamps
            self=varargin{1};
            self.get_MWorks_bitCodes()
            M=self.bitCodes.MWorks_bitCodes;
            sel=M(:,1)>0;
            M=M(sel,:);
            
            
            T=M(:,1);
            T=T-min(T);
            D=diff(T);
            breaks=cat(1,find(D>2.5),length(T));
            
            break_matrix=[breaks(1:end-1) breaks(2:end) diff(breaks)+1];
            break_matrix(break_matrix(:,3)<200,:)=[];
            break_matrix=cat(2,(1:size(break_matrix,1))',break_matrix);
            
            nSessions=size(break_matrix,1);
            
            for iSession=1%1:nSessions
                row=break_matrix(iSession,:);
                timestamps=M(row(2):row(3),1);
                bit_codes=M(row(2):row(3),2);
                
                nFrames=self.mov_info.nFrames;
                frame_rate=self.mov_info.frame_rate;
                
                if 0
                    N=500;
                    T=M(1:N,1);T=T-T(2);
                    [(1:N)' T diff(M(1:N+1,1)) M(1:N,2)]
                    die
                end
                
                output=[timestamps-timestamps(2)+1/frame_rate bit_codes];
                % fix first frame: bitCode is correct, duration not
                output(1,1)=0;
                
                % fix last frame: only runs until nFrames
                output(end,1)=(nFrames-1)/frame_rate;
                
                frame_vector=round(output(:,1)*frame_rate)+1;
                output=[frame_vector output [diff(output(:,1));0]];
                output(:,1:4)
                
                %output(:,4)
                
                die
                
                frame_durations=round(diff(timestamps)*self.mov_info.frame_rate);
                frame_vector=[1 ; [0 ; cumsum(frame_durations)]+2 ; nFrames];
                frame_times=(frame_vector-1)/frame_rate;
                %[frame_times bit_codes [diff(frame_times) ; 0]]
                
                
                
                frame_matrix=zeros(nFrames,8)-1;
                
                
                
                cur_trial=0;
                cur_idx=0;
                timestamps=timestamps-timestamps(1);
                for iFrame=1:nFrames
                    frame_time=(iFrame-1)/frame_rate;
                    if iFrame==1 % No stim in first frame
                        frame_matrix(iFrame,1:2)=[iFrame frame_time];
                    else
                        idx=find(timestamps>frame_time,1,'first');
                        if ~isempty(idx)
                            if cur_idx~=idx
                                %timestamps(idx)
                                cur_idx=idx;
                                cur_trial=cur_trial+1;
                                bit_code=bit_codes(idx);
                                %%% Build scim_bitCodes up from here
                            end
                            frame_matrix(iFrame,1:4)=[iFrame frame_time cur_trial bit_code];
                            
                            %%% Last distractor runs till last frame, next
                            %%% change only when new session starts. so we
                            %%% have to crop it to fill up the nFrames
                        else
                            frame_matrix(iFrame,1:4)=[iFrame frame_time cur_trial+1 M(row(3)+2,2)];
                        end
                    end
                end
                
            end
            
            
        end
        
        
        %%% Get FOV info
        function get_FOV_info(varargin)
            self=varargin{1};
            if nargin>=2
                self.FOV_info.pixel_size_micron=varargin{2};
            else
                self.FOV_info.pixel_size_micron=[.85 .85];
            end
            if self.is_static_FOV()
                if isempty(self.frame_info(1).xyz_submicron)
                    % not there for older files
                    self.FOV_info.coords=self.frame_info(1).xyz_micron;
                else
                    self.FOV_info.coords=self.frame_info(1).xyz_submicron;
                end
                
                self.FOV_info.center=self.FOV_info.coords(1:2);
                self.FOV_info.Z_depth=self.FOV_info.coords(3);
                self.FOV_info.size_px=[self.mov_info.Width self.mov_info.Height];
                self.FOV_info.pixel_aspect_ratio=self.FOV_info.pixel_size_micron(1)/self.FOV_info.pixel_size_micron(2);
                self.FOV_info.size_um=self.FOV_info.size_px.*fliplr(self.FOV_info.pixel_size_micron);
                self.updated=1;
            else
                self.FOV_info.size_px=[self.mov_info.Width self.mov_info.Height];
                self.FOV_info.pixel_aspect_ratio=self.FOV_info.pixel_size_micron(1)/self.FOV_info.pixel_size_micron(2);
                self.FOV_info.size_um=self.FOV_info.size_px.*fliplr(self.FOV_info.pixel_size_micron);
                self.updated=1;
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
                %fprintf('Some axes show more variation (std XYZ=[%3.2f %3.2f %3.2f]) than is acceptable (TH=%3.2f)...\n',[std(M) TH])
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
                
                self.elapsed=toc;
                self.last_action='get_scim_bitCodes';
                self.updated=1;
            else
                disp('Using existing scim_bitCodes...')
            end
        end
        
        function get_MWorks_bitCodes(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.bitCodes.MWorks_bitCodes)
                if nargin>=2
                    F=varargin{2};
                else % if none specified, use first mwk file found in data_folder
                    files=scandir(self.folder_info.data_folder,'.mwk');
                    if isempty(files)
                        error('No .mwk file found in folder')
                    else
                        F=files(1).name;
                    end
                end
                mwk_file_name=fullfile(self.folder_info.data_folder,F);
                
                if ismac
                    %disp('Reading MWK file...')
                    A=getCodecs(mwk_file_name);
                    event_codec=A.codec;
                else
                    mat_file_name=strrep(mwk_file_name,'.mwk','.mat');
                    if exist(mat_file_name,'file')==2
                        load(mat_file_name,'codecs','event_code_strings','event_code_selection','events')
                    else
                        error('No converted .mwk file found for this session...')
                    end
                    
                    event_codec=codecs.codec;
                end
                
                %%% Get stim update events
                tag_name='#stimDisplayUpdate';
                [MW_events,nEvents]=get_events_by_name(mwk_file_name,tag_name);
                
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
                
                self.elapsed=toc;
                self.last_action='get_MWorks_bitCodes';
                self.updated=1;
            else
                disp('Using existing MWorks_bitCodes...')
            end
        end
        
        function find_offset(varargin)
            tic
            self=varargin{1};
            if isempty(self.bitCodes.MWorks_bitCodes)
                error('Bitcodes have to be extracted first using the get_MWorks_bitCodes method...')
            else
                if isempty(self.bitCodes.offset)
                    A=self.bitCodes.scim_bitCodes(:,2);
                    B=self.bitCodes.MWorks_bitCodes(:,2);
                    CC=normxcorr2(A,B);
                    [self.bitCodes.max_val,loc]=max(CC);
                    if self.bitCodes.max_val>.99
                        self.bitCodes.offset=loc-length(A)+1;
                    else
                        T_scim=cat(1,self.frame_info(:).timestamp);
                        offset_temp=loc-length(A)+1
                        length(self.bitCodes.MWorks_bitCodes)
                        T_MWorks=self.bitCodes.MWorks_bitCodes(offset_temp:offset_temp+length(A)-1,1);
                        
                        diff(T_scim)
                        diff(T_MWorks)
                        
                        
                        D_scim=diff(T_scim);
                        subplot(211)
                        %plot((1:length(D_scim))/6.2,D_scim)
                        plot(D_scim)
                        subplot(212)
                        plot(diff(T_MWorks))
                        %[(1:length(A))' B(offset_temp:offset_temp+length(A)-1) A]
                        error('no good match found')
                        self.bitCodes.offset=offset_temp;
                    end
                    
                    self.elapsed=toc;
                    self.last_action='find_offset';
                    self.updated=1;
                else
                    disp('Using existing offset...')
                end
            end
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
                %flyback=double(imread(self.file_name,idx(iFrame),'info',info,'PixelRegion',{[self.mov_info.Height self.mov_info.Height],[1 self.mov_info.Width]}))
                if apply_motion_correction==1
                    offset_ij=self.motion_correction.shift_matrix(iFrame,[4 5]);
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
        
        function g_frames=get_frames_GPU(varargin)
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
            
            try
                TF=existsOnGPU(g_frames); % 0 if deleted, 1 if exists
            catch
                TF=0; % never existed
                disp('g_frames not found on GPU')
            end
            
            if TF==1
                nFrames=size(g_frames,3);
                fprintf('%d frames are already loaded to GPU...\n',nFrames)
            else
                N=length(idx);
                info=imfinfo(self.file_name); % never save info, is huge
                g_frames=zeros([self.mov_info.Height-1 self.mov_info.Width N],'gpuArray');
                for iFrame=1:N
                    frame=double(imread(self.file_name,idx(iFrame),'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                    if apply_motion_correction==1
                        offset_ij=self.motion_correction.shift_matrix(iFrame,[4 5]);
                        frame=offsetIm(frame,offset_ij(1),offset_ij(2),0);
                    else
                        % do nothing
                    end
                    if size(frame,3)==3
                        frame=rgb2gray(frame/256);
                        %frame=frame(:,:,1);
                    end
                    g_frames(:,:,iFrame)=frame;
                end
            end
        end
        
        function find_blank_frames(varargin)
            %%% Check for unusually dark frames, laser power not turned up
            %%% yet.
            tic
            self=varargin{1};
            if isempty(self.mov_info.mean_lum)
                %info=imfinfo(self.file_name); % never save info, is huge
                %self.mov_info.mean_lum=zeros(self.mov_info.nFrames,1);
                %for iFrame=1:self.mov_info.nFrames
                %    frame=double(imread(self.file_name,iFrame,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                %    self.mov_info.mean_lum(iFrame,1)=mean(frame(:));
                %end
                frames=self.get_frames();
                self.mov_info.mean_lum=squeeze(mean(mean(frames,1),2));
                self.mov_info.blank_frames=false(size(self.mov_info.mean_lum));
                
                % Search for longer periods that might cause previous condition to fail
                % low std, low mean
                z_scores=(self.mov_info.mean_lum(1:round(end/4))-mean(self.mov_info.mean_lum(round(end/4):end)))/std(self.mov_info.mean_lum(round(end/4):end));
                self.mov_info.blank_frames(1:find(z_scores<-4,1,'last')+1)=true;
                
                self.elapsed=toc;
                self.last_action='find_blank_frames';
                self.updated=1;
            else
                disp('Using existing blank_frames...')
            end
        end
        
        
        %%% Get experiment info from MWorks
        function get_exp_type(varargin)
            self=varargin{1};
            
            exp_type=[];
            exp_name='';
            
            %%% Try first to read from the MWorks variable exptype directly
            tag_name='ExpType';
            expType_events=get_events_by_name(self.bitCodes.mwk_file_name,tag_name,self.bitCodes.event_codec);
            if ~isempty(expType_events)
                exp_type=mode(cat(1,expType_events.data));
                
                tag_name='ExpName_short';
                expName_events=get_events_by_name(self.bitCodes.mwk_file_name,tag_name,self.bitCodes.event_codec);
                if ~isempty(expType_events)
                    exp_name=expName_events(1).data;
                else
                    exp_name_vector={'RSVP','Retinomapping'};
                    exp_name=exp_name_vector{exp_type};
                end
            else
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
            end
            
            self.Experiment_info.exp_type=exp_type;
            self.Experiment_info.exp_name=exp_name;
            self.updated=1;
        end
        
        function get_MWorks_stimulus_info(varargin)
            self=varargin{1};
            
            if mean(eq(self.bitCodes.scim_bitCodes(:,2),self.bitCodes.MWorks_bitCodes(self.bitCodes.offset:self.bitCodes.offset-1+size(self.bitCodes.scim_bitCodes,1),2)))>.99
                stim_times=self.bitCodes.MWorks_bitCodes(self.bitCodes.offset:self.bitCodes.offset-1+size(self.bitCodes.scim_bitCodes,1),1);
                %stim_times([1 end]); % these time are sufficient to capture all events for this experiment
                
                if ~isempty(self.Experiment_info.exp_type)
                    switch self.Experiment_info.exp_type
                        case 1
                            tag_name='DistractorPresentation_Time';
                            stimDuration_events=get_events_by_name(self.bitCodes.mwk_file_name,tag_name,self.bitCodes.event_codec);
                            stimDuration(1)=stimDuration_events(1).data;
                            tag_name='DistractorIdleTime';
                            blankDuration_events=get_events_by_name(self.bitCodes.mwk_file_name,tag_name,self.bitCodes.event_codec);
                            stimDuration(2)=blankDuration_events(1).data;
                            if isfield(self.Experiment_info,'stimDuration')
                                self.Experiment_info=rmfield(self.Experiment_info,'stimDuration');
                            end
                            self.Experiment_info.stim_duration=stimDuration;
                            
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
                            self.updated=1;
                            
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
            B=mode(reshape(self.bitCodes.scim_bitCodes_raw,self.bitCodes.nBitCodes,[]))'; % BV20150816:this line had a 50 hardcoded as nBitCodes, made sessions with different frame_rate return randomized results
            %B=clean_up_bitCodes_raw();
            
            trial_mapping=self.expand_trial_numbers(B);
            stim_matrix=zeros(self.mov_info.nFrames,10)-1;
            for iTrial=1:length(stimulus_data)
                sel=trial_mapping==iTrial;
                stim_props=[iTrial stimulus_data(iTrial).stim_present stimulus_data(iTrial).condition_nr stimulus_data(iTrial).stim_id stimulus_data(iTrial).position_nr stimulus_data(iTrial).pos_x stimulus_data(iTrial).pos_y stimulus_data(iTrial).size_x stimulus_data(iTrial).size_y stimulus_data(iTrial).rotation];
                stim_matrix(sel,1:length(stim_props))=repmat(stim_props,sum(sel),1);
            end
            stim_matrix=[(1:self.mov_info.nFrames)' stim_matrix];
            self.Experiment_info.stim_matrix=stim_matrix;
            
            self.updated=1;
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
            %fprintf('Completed: %03d%%\n',0)
            fprintf('Running motion detection... ')
            for iFrame=1:self.mov_info.nFrames-1
                cur_frame=double(imread(self.file_name,iFrame,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                next_frame=double(imread(self.file_name,iFrame+1,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                
                diff_vector(iFrame+1)=sum(sqrt((cur_frame(:)-next_frame(:)).^2));
                if ismac
                    del_str=repmat('\b',1,5);
                    fprintf([del_str '%03d%%\n'],round(iFrame/(self.mov_info.nFrames-1)*100))
                end
            end
            fprintf('    Done!\n')
            self.motion_proxy=diff_vector;
            
            self.elapsed=toc;
            self.last_action='do_motion_correction';
        end
        
        function set_smoothing_kernel(varargin)
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                kernel_size=varargin{2};
            else
                kernel_size=[6 6];
            end
            if nargin>=3&&~isempty(varargin{3})
                sigma=varargin{3};
            else
                sigma=.65;
            end
            
            self.motion_correction.kernel=bellCurve2(1,kernel_size/2,[sigma sigma],kernel_size,0);
        end
        
        function set_smoothing_kernel_GPU(varargin)
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                kernel_size=varargin{2};
            else
                kernel_size=[6 6];
            end
            if nargin>=3&&~isempty(varargin{3})
                sigma=varargin{3};
            else
                sigma=.65;
            end
            
            self.motion_correction.kernel=gpuArray(bellCurve2(1,kernel_size/2,[sigma sigma],kernel_size,0));
        end
        
        function find_reference_image(varargin)
            % Implementing motion correction method by Poort et al. 2015 Neuron
            tic
            self=varargin{1};
            
            nBlocks=30;
            block_size=10; % less frames, but more time given our lower sampling rate: 3 vs 30Hz
            nSamples=100;
            
            if isempty(self.motion_correction.reference_image.im)
                fprintf('Finding best reference image...')
                
                %%% Get frames
                frames=self.get_frames();
                
                start_frame=find(self.mov_info.blank_frames==0,1,'first'); % find first non-blank frame
                end_frame=self.mov_info.nFrames;
                
                %%% space 30 blocks of 30 linearly across the stack (or 1s)
                A=floor(linspace(start_frame,end_frame-block_size,nBlocks));
                M=[A(:) A(:)+block_size];
                
                %%% collect candidate reference images
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
                shift_matrix=zeros(nSamples,6,nBlocks);
                for iBlock=1:nBlocks
                    ref=reference_image_candidates(:,:,iBlock);
                    for iSample=1:nSamples
                        sample=sample_images(:,:,iSample);
                        % smooth B?
                        sample=convn(sample,self.motion_correction.kernel,'same');
                        
                        [r,c]=PCdemo(sample,ref);
                        [CC_max,offset]=im_align(sample,ref);
                        
                        shift_matrix(iSample,:,iBlock)=[iSample c r offset CC_max];
                    end
                end
                
                % find best
                total_shift=zeros(nBlocks,1);
                for iBlock=1:nBlocks
                    %total_shift(iBlock)=sum(sqrt(sum(shift_matrix(:,2:3,iBlock).^2,2)));
                    total_shift(iBlock)=sum(sqrt(sum(shift_matrix(:,4:5,iBlock).^2,2)));
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
                self.updated=1;
                
                fprintf(' took %3.2f seconds.\n',self.elapsed)
            else
                disp('Using existing reference image...')
            end
        end
        
        function find_reference_image_GPU(varargin)
            % Implementing motion correction method by Poort et al. 2015 Neuron
            tic
            self=varargin{1};
            
            nBlocks=30;
            block_size=10; % less frames, but more time given our lower sampling rate: 3 vs 30Hz
            nSamples=100;
            
            if isempty(self.motion_correction.reference_image.im)
                fprintf('Finding best reference image...')
                
                %%% Get frames
                %g_frames=gpuArray(self.get_frames(1:30));
                g_frames=self.get_frames_GPU(1);
                
                start_frame=find(self.mov_info.blank_frames==0,1,'first') % find first non-blank frame
                end_frame=self.mov_info.nFrames
                
                %%% space 30 blocks of 30 linearly across the stack (or 1s)
                A=floor(linspace(start_frame,end_frame-block_size,nBlocks));
                M=[A(:) A(:)+block_size];
                
                %%% collect candidate reference images
                g_reference_image_candidates=zeros([size(g_frames,1) size(g_frames,2) nBlocks],'gpuArray');
                for iBlock=1:nBlocks
                    idx=M(iBlock,1):M(iBlock,2);
                    %g_reference_image_candidates(:,:,iBlock)=mean(g_frames(:,:,idx),3);
                    g_reference_image_candidates(:,:,iBlock)=mean(self.get_frames_GPU(idx),3);
                end
                
                die
                % space 100 single frames linearly across the stack
                A=floor(linspace(start_frame,end_frame,nSamples));
                g_sample_images=g_frames(:,:,A);
                
                % cross-correlate
                %self.motion_correction.reference_image.CC_matrix=zeros(nSamples,nBlocks);
                shift_matrix=zeros(nSamples,6,nBlocks);
                for iBlock=1:nBlocks
                    g_ref=g_reference_image_candidates(:,:,iBlock);
                    for iSample=1:nSamples
                        g_sample=g_sample_images(:,:,iSample);
                        % smooth B?
                        g_sample=convn(g_sample,self.motion_correction.kernel,'same');
                        
                        [r,c]=PCdemo(g_sample,g_ref);
                        [CC_max,offset]=im_align(g_sample,g_ref);
                        
                        shift_matrix(iSample,:,iBlock)=[iSample gather(c) gather(r) offset CC_max];
                    end
                end
                
                % find best
                total_shift=zeros(nBlocks,1);
                for iBlock=1:nBlocks
                    %total_shift(iBlock)=sum(sqrt(sum(shift_matrix(:,2:3,iBlock).^2,2)));
                    total_shift(iBlock)=sum(sqrt(sum(shift_matrix(:,4:5,iBlock).^2,2)));
                end
                [min_val,iBest]=min(total_shift);
                
                self.motion_correction.reference_image.idx=M(iBest,1):M(iBest,2);
                self.motion_correction.reference_image.im=g_reference_image_candidates(:,:,iBest);
                self.motion_correction.reference_image.shift_matrix=shift_matrix;
                self.motion_correction.reference_image.min_val=min_val;
                self.motion_correction.reference_image.iBest=iBest;
                self.motion_correction.reference_image.total_shift=total_shift;
                
                self.elapsed=toc;
                self.last_action='find_reference_image';
                self.updated=1;
                
                fprintf(' took %3.2f seconds.\n',self.elapsed)
            else
                disp('Using existing reference image...')
            end
        end
        
        function reset_reference_image(varargin)
            self=varargin{1};
            self.motion_correction.reference_image.im=[];
        end
        
        function do_motion_correction(varargin)
            tic
            self=varargin{1};
            
            if isempty(self.motion_correction.shift_matrix)
                info=imfinfo(self.file_name); % never save info, is huge
                
                ref=self.motion_correction.reference_image.im;
                
                fprintf('Running motion correction... ')
                if ismac
                    str=sprintf('%03d%%\n',0);
                    fprintf(str)
                end
                shift_matrix=zeros(self.mov_info.nFrames,6);
                frames=self.get_frames();
                for iFrame=1:self.mov_info.nFrames
                    %cur_frame=double(imread(self.file_name,iFrame,'info',info,'PixelRegion',{[1 self.mov_info.Height-1],[1 self.mov_info.Width]}));
                    cur_frame=frames(:,:,iFrame);
                    if size(cur_frame,3)==3
                        cur_frame=rgb2gray(cur_frame/256);
                    end
                    
                    %%% smooth frame
                    kernel_size=[6 6];
                    sigma=.65;
                    kernel=bellCurve2(1,kernel_size/2,[sigma sigma],kernel_size,0);
                    cur_frame=convn(cur_frame,kernel,'same');
                    
                    %%% Calc normxcorr2
                    [CC_max,offset]=im_align(cur_frame,ref);
                    
                    %%% Calc phase shift
                    %[r,c]=PCdemo(cur_frame,ref);
                    r=-1;
                    c=-1;
                    %[iFrame c r]
                    shift_matrix(iFrame,:)=[iFrame c r offset CC_max];
                    
                    
                    %%% Raise error when mismatch
                    if ismac
                        del_str=repmat('\b',1,length(str)-2);
                        fprintf([del_str '%03d%%\n'],round(iFrame/(self.mov_info.nFrames-1)*100))
                    end
                end
                fprintf('Done! (Took %3.2f seconds)\n',toc)
                self.motion_correction.shift_matrix=shift_matrix;
                
                self.elapsed=toc;
                self.last_action='do_motion_correction';
                self.updated=1;
            else
                disp('Using existing shift_matrix...')
            end
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
            
            M=self.motion_correction.shift_matrix;
            d=calc_dist([M(:,4:5) repmat(mean(M(:,4:5)),size(M,1),1)]); % euclidean distance moved vs center point, in px
            self.motion_correction.motion_frames=d>5; % too large shift in FOV, cells might get out of the frame
            self.motion_correction.mismatch_frames=M(:,6)<.65; % max crosscorrelation too low, probably z-motion
            
            %T=self.motion_correction.shift_matrix(:,1);
            D=calc_dist([self.motion_correction.shift_matrix(:,4:5)*0 self.motion_correction.shift_matrix(:,4:5)]);
            M=[cat(1,D,[0; 0;0]) cat(1,0,D,[0;0]),cat(1,[0; 0],D,0),cat(1,[0; 0;0],D)];M(end-2:end,:)=[];
            self.motion_correction.ignore_frames=std(M,[],2)>TH;
            self.motion_correction.variability_threshold=TH; % too much variability, probably heavy motion artifacts
            
            self.elapsed=toc;
            self.last_action='find_motion_frames';
            self.updated=1;
        end
        
        
        %%% Calc MIPs using motion correction
        function do_calc_MIPs(varargin)
            tic
            self=varargin{1};
            
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
                    if ~ismac % only if we are on the server, takes a long time to run
                        self.MIP_cc_local.data=CrossCorrImage(frames);
                        self.MIP_cc_local.gamma_val=.6;
                    end
                end
                if isempty(self.MIP_kurtosis.data)
                    self.MIP_kurtosis.data=3-kurtosis(frames,[],3);
                    self.MIP_kurtosis.gamma_val=.6;
                end
                fprintf('Done!\n')
                
                self.elapsed=toc;
                self.last_action='do_calc_MIPs';
                self.updated=1;
            else
                disp('Using existing MIPs...')
            end
        end
        
        function reset_MIPs(varargin)
            self=varargin{1};
            self.MIP_avg.data=[];
            self.MIP_max.data=[];
            self.MIP_std.data=[];
        end
        
        
        %%% Find ROIs
        function find_ROIs(varargin)
            self=varargin{1};
            
            tic
            
            window_size=[40 40];
            
            %%% Load MIP
            IM=self.MIP_cc_local.data;
            
            %%% Adaptive threshold
            try
                ws=20;
                C=-.03;
                tm=0;
                bw=adaptivethreshold(IM,ws,C,tm);
            catch
                self.imshow(IM)
            end
            
            %%% Remove speckle
            SE=[0 1 0 ; 1 1 1 ; 0 1 0];
            bw=imopen(bw,SE);
            
            %%% Filter based on spot size
            props=regionprops(bw,'Area','Centroid');
            coords=cat(1,props.Centroid);
            TH.area=[50 500];
            sel_area=between(cat(1,props.Area),TH.area);
            
            %%% Filter based on distance from the edge
            sel_edge=all([between(coords(:,2),[window_size(1)/2+1 self.mov_info.Height-window_size(1)/2]) between(coords(:,1),[window_size(2)/2+1 self.mov_info.Width-window_size(2)/2])]')';
            
            coords=coords(sel_area&sel_edge,:);
            bw=bwselect(bw,coords(:,1),coords(:,2));
            
            %%% Convert into ROIs
            bw_separated=bwlabel(bw);
            ROI_vector=unique(bw_separated(:));
            ROI_vector=ROI_vector(ROI_vector>0);
            nROI=length(ROI_vector);
            
            for iROI=1:nROI
                ROI_nr=ROI_vector(iROI);
                sel=bw_separated==ROI_nr;
                neg=imerode(sel,SE);
                edge=sel-neg;
                
                [X,Y]=find(edge);
                
                ellipse_properties=fit_ellipse(Y,X);
                
                blank=self.blank_ROI();
                ROI=blank.ROI;
                
                ROI.ROI_nr=ROI_nr;
                ROI.base_coord=[ellipse_properties.X0_in ellipse_properties.Y0_in];
                ROI.nCoords=size(coords,1);
                ROI.coords=[Y-ROI.base_coord(1)+window_size(1)/2+1 X-ROI.base_coord(2)+window_size(2)/2+1];
                
                ROI.ellipse_properties=ellipse_properties;
                ROI.ellipse_coords=[ellipse_properties.rotated_ellipse(1,:)'-ROI.base_coord(1)+window_size(1)/2+1 ellipse_properties.rotated_ellipse(2,:)'-ROI.base_coord(2)+window_size(2)/2+1];
                ROI.ellipse_coords_centered=ROI.ellipse_coords;
                ROI.coords_MIP=ellipse_properties.rotated_ellipse';
                ROI.coords_MIP_plot=ellipse_properties.rotated_ellipse';
                
                ROI.center_coords=[ellipse_properties.X0_in ellipse_properties.Y0_in];
                
                ROI.ROI_rect=round([ROI.center_coords ROI.center_coords]+[-window_size/2+1 window_size/2]);
                
                %%% Create masks
                % Get region from image around ROI
                try
                    T=IM(ROI.ROI_rect(2):ROI.ROI_rect(4),ROI.ROI_rect(1):ROI.ROI_rect(3));
                catch
                    figure()
                    self.imshow(sel)
                    title('Failed to fit ROI')
                    break
                end
                
                % Generate mask to isolate soma pixels
                %offset=ROI.ellipse_coords_centered(1,:)-window_size/2
                mask_soma=poly2mask(ROI.ellipse_coords_centered(:,1),ROI.ellipse_coords_centered(:,2),size(T,1),size(T,2));
                %figure(50)
                %ROI.ellipse_coords_centered
                %imshow(mask_soma,[])
                
                %die
                
                % Generate mask to isolate neuropil pixels
                mask_neuropil=imresize(mask_soma,sqrt(2),'nearest');
                ext=round((size(mask_neuropil)-window_size)/2);
                mask_neuropil=mask_neuropil(ext:ext+window_size(1)-1,ext:ext+window_size(2)-1);
                mask_neuropil=mask_neuropil-mask_soma;
                
                ROI.mask_soma=mask_soma;
                ROI.mask_neuropil=mask_neuropil;
                
                ROI.timeseries_soma=[];
                %ROI.time_series_neuropil=[];
                
                %if isfield(self.ROI_definitions(1).ROI,'timeseries_neuropil')
                %    self.ROI_definitions(1).ROI=rmfield(self.ROI_definitions(1).ROI,'timeseries_neuropil');
                %end
                %%% Store ROI properties
                try
                    self.ROI_definitions(1).ROI(iROI)=ROI;
                catch
                    self.ROI_definitions(1).ROI(iROI)
                    ROI
                    error('Objects are not identical...')
                end
            end
            fprintf('Found %d ROIs!\n',nROI)
            
            self.elapsed=toc;
            self.last_action='find_ROIs';
            self.updated=1;
        end
        
        function reset_ROIs(varargin)
            self=varargin{1};
            self.ROI_definitions(1)=self.blank_ROI();
        end
        
        function blank_ROI=blank_ROI(varargin)
            blank_ROI=struct('ROI',struct(...
                'ROI_nr',[],'base_coord',[],'nCoords',[],'coords',[],...
                'ellipse_properties',[],'ellipse_coords',[],'coords_MIP',[],'coords_MIP_plot',[],'center_coords',[],'ellipse_coords_centered',[],...
                'ROI_rect',[],'mask_soma',[],'mask_neuropil',[],'timeseries_soma',[])...
                );
            
        end
        
        function plot_ROIs(varargin)
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                ROI_definition_nr_local=varargin{2};
            else
                ROI_definition_nr_local=self.ROI_definition_nr;
            end
            im=self.MIP_cc_local;
            ROI=self.ROI_definitions(ROI_definition_nr_local).ROI;
            nROIs=length(ROI);
            
            self.imshow(im)
            hold on
            for iROI=1:nROIs
                coords=ROI(iROI).coords_MIP;
                plot(coords(:,1),coords(:,2),'r')
            end
            hold off
        end
        
        
        %%% Do ROI extraction
        function create_mask_from_ROI(varargin)
            self=varargin{1};
            
            ROIs=get_ROI_definitions(self);
            N=length(ROIs);
            
            %%% loop tru ROIs and create sparse matrix with masks
            for iROI=1:N
                FOV_rect_px=self.FOV_info.size_px;
                
                stretch_factor=round(self.FOV_info.pixel_size_micron/min(self.FOV_info.pixel_size_micron));
                stretch_coords=round(fliplr(self.FOV_info.size_px).*stretch_factor);
                
                ROI=ROIs(iROI);
                
                %%% Construct mask for soma
                mask_soma=poly2mask(ROI.coords_MIP(:,1),ROI.coords_MIP(:,2),self.FOV_info.size_px(2),self.FOV_info.size_px(1));
                
                %%% Construct mask for neuropil
                mask_soma=imresize(mask_soma,stretch_coords,'nearest');
                mask_neuropil=bwmorph(mask_soma,'thicken',4)-mask_soma;
                
                %%% Resize to original aspect ratio
                mask_soma=imresize(mask_soma,fliplr(FOV_rect_px),'nearest');
                mask_neuropil=imresize(mask_neuropil,fliplr(FOV_rect_px),'nearest');
                
                %%% Store data
                self.ROI_matrix(iROI).mask_soma=sparse(mask_soma);
                self.ROI_matrix(iROI).mask_neuropil=sparse(mask_neuropil);
                self.ROI_matrix(iROI).data=[iROI ROI.ROI_nr sum(mask_neuropil(:)) sum(mask_soma(:)) sum(mask_neuropil(:))/sum(mask_soma(:))];
                
                %%% Add bounding box info
                P=regionprops(mask_soma,'boundingBox');
                self.ROI_matrix(iROI).rect_soma=floor(P.BoundingBox);
                P=regionprops(mask_neuropil,'boundingBox');
                self.ROI_matrix(iROI).rect_neuropil=floor(P.BoundingBox);
            end
            
            if 0
                %% Vizualize using this code
                mask_soma=combine_sparse_masks(self.ROI_matrix,'mask_soma');
                mask_neuropil=combine_sparse_masks(self.ROI_matrix,'mask_neuropil');
                
                figure(132)
                imshow(imresize(mask_neuropil*.5+mask_soma,stretch_coords,'nearest'),[])
            end
        end
        
        function clean_neuropil_shell(varargin)
            % This method removes all pixels from the neuropil shell that
            % overlap with neighboring cells
            self=varargin{1};
            ROIs=get_ROI_definitions(self);
            N=length(ROIs);
            mask_soma=combine_sparse_masks(self.ROI_matrix,'mask_soma');
            for iROI=1:N
                hole=(self.ROI_matrix(iROI).mask_neuropil+mask_soma) > 1; % find overlapping pixels
                self.ROI_matrix(iROI).mask_neuropil(hole==1)=0; % reject overlapping pixels
            end
        end
        
        function do_trace_extraction(varargin)
            tic
            self=varargin{1};
            
            temp=get_ROI_definitions(self);
            if isempty(self.Activity_traces.activity_matrix)||length(temp)~=size(self.Activity_traces.activity_matrix,2)
                %% Grab average ROI activity over the motion corrected movie
                %if ~isfield(self.ROI_definitions(ROI_definition_nr).ROI,'timeseries_soma')
                
                if isempty(temp)
                    disp('No valid ROI definitions found...')
                else
                    if ~isfield(temp(end),'timeseries_soma')||isempty(temp(end).timeseries_soma)
                        ROIs=self.extract_traces();
                        self.ROI_definitions(self.ROI_definition_nr).ROI=ROIs;
                    end
                end
            else
                % load existing definitions
                ROIs=get_ROI_definitions(self);
                disp('Using existing trace_matrix')
            end
            
            %ROI_vector=cat(1,ROIs.ROI_nr);
            nROI=length(ROIs);
            nFrames=self.mov_info.nFrames;
            
            activity_matrix=zeros(nROI,nFrames);
            spike_matrix=activity_matrix;
            
            fprintf('Using extraction method #%d\n',self.Activity_traces.extraction_options.calc_delta_f_method)
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
            self.updated=1;
            
            
            
        end
        
        function reset_trace_matrix(varargin)
            self=varargin{1};
            self.Activity_traces.activity_matrix=[];
            temp=self.ROI_definitions(self.ROI_definition_nr).ROI;
            nROI=length(temp);
            for iROI=1:nROI
                self.ROI_definitions(self.ROI_definition_nr).ROI(iROI).timeseries_soma=[];
                self.ROI_definitions(self.ROI_definition_nr).ROI(iROI).timeseries_neuropil=[];
            end
        end
        
        function ROIs=extract_traces(varargin)
            %%% Goes way faster than before!
            self=varargin{1};
            
            tic
            fprintf('Loading frames... ')
            %self.file_name
            frames=self.get_frames([],1); % get all frames, apply motion correction
            %frames=self.get_frames(1:20,1); % get all frames, apply motion correction
            F=frames(:);
            H=size(frames,1);
            W=size(frames,2);
            nFrames=size(frames,3);
            
            %nFrames=size(frames,3);
            fprintf('took: %3.1fs\n',toc);
            
            ROIs=get_ROI_definitions(self);
            nROI=length(ROIs);
            tic
            fprintf('Extracting %d ROIs data... ',nROI)
            for iROI=1:nROI
                rect=ROIs(iROI).ROI_rect;
                
                %%% Allow selection of ROI close to border, pad with zeros
                %%% if over!
                switch 3
                    case 1
                        tic
                        vol=frames(rect(2):rect(4),rect(1):rect(3),:);
                        
                        %%% Extract soma pixels
                        mask=repmat(ROIs(iROI).mask_soma,1,1,size(vol,3));
                        res=vol.*mask;
                        ROIs(iROI).timeseries_soma=squeeze(mean(mean(res,1),2));
                        
                        %%% Extract neuropil pixels
                        mask=repmat(ROIs(iROI).mask_neuropil,1,1,size(vol,3));
                        res=vol.*mask;
                        ROIs(iROI).timeseries_neuropil=squeeze(mean(mean(res,1),2));
                        %squeeze(mean(mean(res,1),2))
                        
                        %figure()
                        %subplot(221)
                        %self.imshow(mean(vol,3))
                        %subplot(222)
                        %self.imshow(mean(res,3))
                        %subplot(223)
                        %self.imshow(mean(res,3))
                        toc
                    case 2
                        tic
                        % create and store outside contour of soma
                        % upon extraction, fill this and get indices in
                        % whole image space
                        % crop indices outside of image
                        % repeat indices for nFrames
                        % extract pixels
                        % reshape and average per image plane
                        
                        %vol=frames(rect(2):rect(4),rect(1):rect(3),:);
                        w=size(ROIs(iROI).mask_soma,1)/2;
                        
                        % get indices in image space
                        [x,y]=find(ROIs(iROI).mask_soma);
                        x=x+round(ROIs(iROI).center_coords(1)-w);
                        y=y+round(ROIs(iROI).center_coords(2)-w);
                        
                        % crop indices
                        x(x<0)=0;x(x>W)=W;
                        y(y<0)=0;y(y>H)=H;
                        
                        % vectorize
                        indices=y*W+x;
                        A=zeros(W*H,1);
                        A(indices)=1;
                        
                        % repeat
                        I=repmat(A,nFrames,1);
                        
                        % extract pixels
                        pixels=F(I==1);
                        
                        % reshape and average
                        response=reshape(pixels,[],nFrames);
                        
                        % write to data
                        ROIs(iROI).timeseries_soma=mean(response)';
                        size(mean(response)')
                        toc
                    case 3
                        % use newly created sparse matrix as input,
                        % create bounding box
                        % extract only that from the 3D frames matrix
                        % now replicated the selection plane to a matching
                        % 3D matrix
                        
                        rect_soma=self.ROI_matrix(iROI).rect_soma;
                        mask_soma=self.ROI_matrix(iROI).mask_soma;
                        rect_neuropil=self.ROI_matrix(iROI).rect_neuropil;
                        mask_neuropil=self.ROI_matrix(iROI).mask_neuropil;
                        
                        %%% Get time series for soma pixels
                        indices_X=rect_soma(2)+1:rect_soma(2)+rect_soma(4)-1;
                        indices_Y=rect_soma(1)+1:rect_soma(1)+rect_soma(3)-1;
                        vol=frames(indices_X,indices_Y,:);
                        mask=repmat(full(mask_soma(indices_X,indices_Y)),1,1,nFrames);
                        res=vol.*mask;
                        ROIs(iROI).timeseries_soma=squeeze(mean(mean(res,1),2));
                        
                        %%% Get time series for soma pixels
                        indices_X=rect_neuropil(2)+1:rect_neuropil(2)+rect_neuropil(4)-1;
                        indices_Y=rect_neuropil(1)+1:rect_neuropil(1)+rect_neuropil(3)-1;
                        vol=frames(indices_X,indices_Y,:);
                        mask=repmat(full(mask_neuropil(indices_X,indices_Y)),1,1,nFrames);
                        res=vol.*mask;
                        ROIs(iROI).timeseries_neuropil=squeeze(mean(mean(res,1),2));
                        
                        if 0
                            %%
                            figure(123)
                            imshow(mean(res,3),[])
                            colormap(self.green)
                        end
                end
                
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
                %[mean(F_raw) mean(F_neuropil)]
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
                case 0
                    delta_F_no_drift=F_no_drift; % leave signal as is
                    y_label='raw F';
                    fixed_y_scale=10000;
                case 1 % naive way
                    delta_F_no_drift=(F_no_drift-mean(F_no_drift))/mean(F_no_drift);
                    y_label='\DeltaF/F';
                    fixed_y_scale=30;
                case 2 % Feinberg
                    delta_F_no_drift=F_no_drift/mean(F_no_drift)-1;
                    y_label='\DeltaF/F';
                    fixed_y_scale=30;
                case 3 % z-score naive
                    delta_F_no_drift=(F_no_drift-mean(F_no_drift))/std(F_no_drift);
                    y_label='Z(F)';
                    fixed_y_scale=30;
                case 4 % minimal std for z-score
                    disp('Under construction')
                    delta_F_no_drift=(F_no_drift-median(F_no_drift));
                    
                    % estimate sigma, minimum from 10s sliding window sigma's
                    averaging_interval=10; % seconds
                    nFrames_avg=averaging_interval*round(frame_rate);
                    if rem(nFrames_avg,2)==1
                        nFrames_avg=nFrames_avg+1;
                    end
                    
                    sigma_vector=zeros(size(delta_F_no_drift));
                    for iSample=1:nFrames
                        if between(iSample,[nFrames_avg/2+1 nFrames-nFrames_avg/2])
                            sigma_values=delta_F_no_drift(iSample-nFrames_avg/2+1:iSample+nFrames_avg/2);
                        elseif iSample<nFrames_avg
                            sigma_values=delta_F_no_drift(1:iSample+nFrames_avg/2);
                        elseif iSample>nFrames-nFrames_avg
                            sigma_values=delta_F_no_drift(iSample-nFrames_avg/2+1:end);
                        else
                            die
                        end
                        sigma_vector(iSample)=std(sigma_values);
                    end
                    sigma_est=min(sigma_vector);
                    delta_F_no_drift=delta_F_no_drift/sigma_est;
                    y_label='Z(F)';
                    fixed_y_scale=30;
                case 5 % z-score based on e-phys analysis
                    % find regions without 'spike' and use rest for std and mean calculation
                    F_temp=F_no_drift-median(F_no_drift);
                    noiseSTD=median(abs(F_temp))/.6745;
                    min_heigth_factor=4;
                    TH=noiseSTD*min_heigth_factor;
                    selection_vector=F_temp<TH;
                    mu=mean(F_no_drift(selection_vector));
                    sigma=std(F_no_drift(selection_vector));
                    
                    % Store normalizaion values
                    self.Activity_traces.normalization_matrix(iROI,1:4)=[iROI extraction_options.calc_delta_f.method mu sigma];
                    
                    delta_F_no_drift=(F_no_drift-mu)/sigma;
                    if 0
                        %%
                        T=1:length(F_temp);
                        plot(T,F_temp)
                        hold on
                        plot([0 length(F_temp)],[TH TH])
                        plot(T(selection_vector==1),F_temp(selection_vector==1),'r')
                        hold off
                    end
                    y_label='Z(F)';
                    fixed_y_scale=30;
                    
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
            self.Activity_traces.extraction_options=extraction_options;
            %self.Activity_traces.normalization_matrix=normalization_matrix;
            
            delta_F_no_drift(blank_frames)=0;
            
            trace=delta_F_no_drift;
        end
        
        
        %%% Join datasets
        function [clusters,cluster_vector,nClusters,files]=find_FOV_clusters(varargin)
            %self=varargin{1};
            folder=varargin{2};
            if nargin>=3&&~isempty(varargin{3})
                min_dist=varargin{3};
            else
                min_dist=20;
            end
            
            files=scandir(folder,'mat');
            nFiles=length(files);
            center_coords=zeros(nFiles,3);
            valid_FOV=zeros(nFiles,1);
            for iFile=1:nFiles
                load_name=fullfile(folder,files(iFile).name);
                load(load_name,'session_data')
                
                if session_data.is_static_FOV()&&session_data.mov_info.nFrames>300
                    center_coords(iFile,:)=session_data.FOV_info.coords;
                    valid_FOV(iFile)=1;
                else
                    valid_FOV(iFile)=0;
                end
            end
            Z=linkage(center_coords,'average');
            C=cluster(Z,'cutoff',min_dist,'Criterion','distance');
            
            % renumber clusters
            cluster_vector=unique(C,'stable');
            clusters=zeros(size(C));
            nClusters=length(cluster_vector);
            for iC=1:nClusters
                sel=C==cluster_vector(iC);
                clusters(sel)=iC;
            end
            clusters(valid_FOV==0)=NaN;
            %cluster_numbers=1:size(C,1);
            cluster_vector=unique(clusters(~isnan(clusters)));
            nClusters=length(unique(clusters(~isnan(clusters))));
        end
        
        function ROI_counts=get_ROI_counts(varargin)
            self=varargin{1};
            nSessions=length(self);
            ROI_counts=zeros(nSessions,1);
            for iSession=1:nSessions
                ROIs=get_ROI_definitions(self(iSession));
                %ROIs=self(iSession).ROI_definitions(self.ROI_definition_nr).ROI;
                ROI_counts(iSession)=length(ROIs);
            end
        end
        
        function distances=get_ROI_distances(varargin)
            self=varargin{1};
            ROI_definition_nr=varargin{2};
            nSessions=length(self);
            if nSessions>1
                ref=self(1).ROI_definitions(ROI_definition_nr).ROI;
                nROIs=length(ref);
                distances=zeros(nROIs,nSessions-1);
                for iSession=2:nSessions
                    ROIs=self(iSession).ROI_definitions(ROI_definition_nr).ROI;
                    A=cat(1,ref.center_coords);
                    B=cat(1,ROIs.center_coords);
                    distances(:,iSession-1)=calc_dist([A B]);
                end
            else
                distances=[];
            end
        end
        
        function dataset=join_data_sessions(varargin)
            if nargin>=1
                self=varargin{1};
                
                dataset=imaging_datasets(self(1),self.ROI_definition_nr);
                nSessions=length(self);
                STIM=[];
                RESP=[];
                SPIKE=[];
                for iSession=1:nSessions
                    M=self(iSession).Experiment_info.stim_matrix;
                    R=self(iSession).Activity_traces.activity_matrix;
                    S=self(iSession).Activity_traces.spike_matrix;
                    
                    %%% Clip last trial or last trial and blank
                    nTrials=max(M(:,2));
                    if mod(nTrials,2)==1
                        nTrials_use=nTrials-2;
                    else
                        nTrials_use=nTrials-1;
                    end
                    sel=M(:,2)<=nTrials_use & M(:,3)>=0;
                    
                    %%% 2DO: add selection based on blank_frames, ignore frames
                    %%% etc...
                    
                    STIM=cat(1,STIM,M(sel,:));
                    RESP=cat(1,RESP,R(sel,:));
                    SPIKE=cat(1,SPIKE,S(sel,:));
                end
                %% reformat cols 1 and 2
                STIM(:,1)=1:size(STIM,1);
                %trial_nrs=STIM(:,2);
                C=STIM(:,3); C(C==0)=-1;
                condition_matrix=parse_conditions(C);
                nTrials=max(condition_matrix(:,1));
                STIM(:,2)=0;
                for iTrial=1:nTrials
                    row=condition_matrix(iTrial,:);
                    if iTrial<nTrials
                        next_row=condition_matrix(iTrial+1,:);
                        STIM(row(2):next_row(3),2)=row(1);
                    else
                        STIM(row(2):size(STIM,1),2)=row(1);
                    end
                end
                
                %%% Fill dataset
                dataset.exp_type=self(1).Experiment_info.exp_type;
                dataset.exp_name=self(1).Experiment_info.exp_name;
                dataset.stim_duration=self(1).Experiment_info.stim_duration;
                dataset.STIM=STIM;
                dataset.RESP=RESP;
                dataset.SPIKE=SPIKE;
                frame_time=1/self(1).mov_info.frame_rate;
                dataset.timeline=dataset.STIM(:,1)*frame_time;
            end
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
            %self.analyse_data_new([Stim(:,1:3) Stim(:,variable) zscore(self.mov_info.mean_lum)],iROI)
            
            %self.analyse_data([Stim(:,1:3) Stim(:,variable) Resp(:,iROI)])
            self.analyse_data_new([Stim(:,1:3) Stim(:,variable) Resp(:,iROI)],iROI)
        end
        
        function analyse_data_new(varargin)
            self=varargin{1};
            data=varargin{2};
            iROI=varargin{3};
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
                [s,mess,messid]=mkdir(self.folder_info.save_folder);
            end
            if self.updated==1
                self.updated=0; % reset updated flag
                session_data=self;
                save(self.save_name,'session_data')
                
                fprintf('Data saved to %s!\n',self.save_name)
                
            else
                disp('No changes detected, not saving...')
            end
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
                if nargin>=3&&~isempty(varargin{3})
                    gamma_val=varargin{3};
                else
                    gamma_val=.5;
                end
            end
            
            if nargin>=4&&~isempty(varargin{4})
                fH=varargin{4};
            else
                fH=444;
            end
            
            nFrames=size(im,3);
            figure(fH)
            %H=imshow(real(calc_gamma(im(:,:,1),gamma_val)),[]);
            H=imagesc(real(calc_gamma(im(:,:,1),gamma_val)));
            pbaspect([self.FOV_info.size_um 1])
            axis off
            colormap(self.green)
            if nFrames>1
                %%% Show movie
                for iFrame=1:nFrames
                    frame=calc_gamma(im(:,:,iFrame),gamma_val);
                    set(H,'Cdata',real(frame))
                    drawnow
                end
            end
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
            p(1)=plot(center(1),center(2),'mo','markersize',23);
            p(2)=plot(center(1),center(2),'ws','markersize',20);
            
            colormap(self.green)
            hold off
            for iFrame=1:self.mov_info.nFrames
                cur_frame=frames(:,:,iFrame);
                coords=center-self.motion_correction.shift_matrix(iFrame,2:3);
                coords2=center-self.motion_correction.shift_matrix(iFrame,4:5);
                
                set(H,'Cdata',calc_gamma(cur_frame,.5))
                set(p(1),'xData',coords(1),'yData',coords(2))
                set(p(2),'xData',coords2(1),'yData',coords2(2))
                if self.motion_correction.ignore_frames(iFrame)==2
                    set(p(1),'marker','o','color','m','markersize',5)
                    set(p(2),'marker','s','color','w','markersize',2)
                else
                    set(p(1),'marker','o','color','m','markersize',23)
                    set(p(2),'marker','s','color','w','markersize',20)
                end
                drawnow
            end
        end
        
        function plot_FOV(varargin)
            self=varargin{1};
            S=self;
            
            nSessions=length(S);
            figure(23)
            clf
            hold on
            for iSession=1:nSessions
                session_data=S(iSession);
                center=session_data.FOV_info.center;
                FOV_rect=[0 0 session_data.FOV_info.size_um];
                ROI=CenterRectOnPoint(FOV_rect,center(1),center(2))/1000;
                %name=session_data.folder_info.raw_name;
                %name=strrep(name,'_',' ');
                
                circle([0 0],2,100,'r-',2);
                plotRect(ROI,'k');
                text(ROI(1)+.05,ROI(2)+.25,sprintf('Depth %3.1fµm',session_data.FOV_info.Z_depth))
                text(ROI(1)+.05,ROI(2)+.1,sprintf('#%d',iSession))
            end
            hold off
            axis square
            axis equal
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
                axis([1 size(A,1) self.Activity_traces.extraction_options.fixed_y_scale*-.1 self.Activity_traces.extraction_options.fixed_y_scale])
                title(sprintf('ROI #%d',iROI))
                set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
            end
        end
        
        function CC=corr_traces(varargin)
            self=varargin{1};
            A=self.Activity_traces.activity_matrix;
            CC=corr(A);
            CC_zero_diag=CC;
            CC_zero_diag(eye(size(CC))==1)=0;
            CC_vector=squareform(CC_zero_diag);
            if usejava('jvm')
                if any(CC_vector>.80)
                    figure(555)
                    plot(sort(CC_vector,'descend'))
                    axis([1 length(CC_vector) -1 1])
                    box off; axis square
                    die
                else
                    figure(555)
                    subplot(121)
                    imagesc(CC)
                    box off; axis square
                    
                    subplot(122)
                    plot(sort(CC_vector,'descend'))
                    line([1 length(CC_vector)],[0 0])
                    axis([1 length(CC_vector) -1 1])
                    box off; axis square
                end
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
        
        function import_movie(varargin)
            self=varargin{1};
            filename=varargin{2};
            
            obj = Tiff(filename,'r');
            
            if obj.isTiled
                obj.numberOfTiles
            else
                obj.numberOfStrips
            end
            obj.PlanarConfiguration
            
            im=double(obj.read);
            self.imshow(im(1:end,:))
        end
        
        function export_movie(varargin)
            % option to do motion correction first
            %self=varargin{1};
            self=varargin{1};
            if nargin>=2&&~isempty(varargin{2})
                filename=varargin{2};
            else
                [f,fn,ext]=fileparts(self.file_name);
                filename=fullfile(f,'Processed',[fn '_motionCorrected' ext]);
                savec(filename)
            end
            if nargin>=3&&~isempty(varargin{3})
                frames=varargin{3};
            else
                frames=self.get_frames([],1); % usually correct motion artefacts
            end
            Tiff(filename,'w');
            
            nFrames=size(frames,3);
            for iFrame=1:nFrames
                obj = Tiff(filename,'a');
                obj.setTag('ImageLength',size(frames,1));
                obj.setTag('ImageWidth', size(frames,2));
                obj.setTag('Photometric', Tiff.Photometric.MinIsBlack);
                obj.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
                obj.setTag('BitsPerSample', 16); % was 8 in example
                obj.write(uint16(frames(:,:,iFrame)));
                obj.close();
            end
        end
        
        function loop_status(varargin)
            self=varargin{1};
            
        end
    end
    
end