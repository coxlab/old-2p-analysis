classdef data_object < handle
    properties
        system=struct;
        dirs=struct;
        files=struct;
        ROI_parameters=struct;
        colormaps=struct;
        
        reconstructed=0;
        registered=0;
        
        save_it=1;
        updated=false;
        log=[];
    end
    
    
    methods
        %%% Constructor
        function self=data_object(varargin)
            t0=clock;
            exp_name=varargin{1};
            
            %%% Get script folder
            self.dirs.script_dir=up1(fileparts(mfilename('fullpath')));
            
            %%% Add subfolder to path
            addpath(genpath(self.dirs.script_dir))
            
            %%% Get info about the environment
            [self.system.host,self.system.user,self.system.platform]=get_ID();
            
            
            switch self.system.host
                case 'dixie'
                    switch self.system.user
                        case 'benvermaercke'
                            root_folder='/nas/volume1/2photon';
                            self.ROI_parameters.ROI_definition_nr=2;
                        case ''
                        otherwise
                            error('Username unknown')
                    end
                case 'BKRUNCH'
                    switch self.system.user
                        case 'ben'
                            root_folder='Q:\data\'; % quadraraid
                            self.ROI_parameters.ROI_definition_nr=2;
                        case ''
                        otherwise
                            error('Username unknown')
                    end
                case 'COXLAB-2P'
                    switch self.system.user
                        case 'labuser'
                            root_folder='D:\Dropbox (coxlab)'; % seagate external drive
                        otherwise
                            error('Username unknown')
                    end
                case 'BEN-PC'
                    switch self.system.user
                        case 'LBP'
                            root_folder='C:\Users\LBP\Dropbox (coxlab)';
                    end
                    
                    
                case 'Bens-MacBook-Pro.local'
                    switch self.system.user
                        case 'benvermaercke'
                            root_folder='/Users/benvermaercke/Dropbox (coxlab)/';
                            self.ROI_parameters.ROI_definition_nr=2;
                        otherwise
                            error('Username unknown')
                    end
                    
                    %case ''
                otherwise
                    error('Hostname unknown...')
            end
            
            %%% Set folder structure
            self.dirs.root_folder=root_folder;
            self.dirs.raw_folder=fullfile(root_folder,'2p-data/raw',exp_name);
            self.dirs.rec_folder=fullfile(root_folder,'2p-data/rec',exp_name);
            self.dirs.reg_folder=fullfile(root_folder,'2p-data/reg',exp_name);
            self.dirs.MWorks_folder=fullfile(root_folder,'MWorks',exp_name);
            self.dirs.eyeTracker_folder=fullfile(root_folder,'eyetracker',exp_name);
            self.dirs.analysis_folder=fullfile(root_folder,'analysis',exp_name);
            
            if isdir(self.dirs.raw_folder)==0
                warning('Data folder not found on this machine...')
            else
                self.files=scandir(self.dirs.raw_folder,'.tif');
                fprintf('Folder found and contains %d files! \n',length(self.files))
            end
            
            %%% Set popular look-up tables
            self.colormaps.red=[linspace(0,1,256)' zeros(256,1) zeros(256,1)];
            self.colormaps.green=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];
            
            self.log=log_object(); % set up log object
            self.updated=true;
            self.log.append('Constructor',t0)
        end
        
        
        function preprocess_data(self,varargin)
            t0=clock;
            nFiles=length(self.files);
            
            warning off
            for iFile=1:nFiles
                tif_name=fullfile(self.dirs.raw_folder,self.files(iFile).name);
                tgt_name=fullfile(self.dirs.rec_folder,self.files(iFile).name);
                MW_name=fullfile(self.dirs.MWorks_folder,'flyback_data.mat');
                savec(tgt_name)
                savec(MW_name)
                if exist(tif_name,'file')==2
                    info=imfinfo(tif_name);
                    nFrames=length(info);
                    
                    if self.save_it==true
                        flybacks=zeros(nFrames,info(1).Width);
                        src=Tiff(tif_name,'r');
                        Tiff(tgt_name,'w');
                        
                        for iFrame=1:nFrames
                            src.setDirectory(iFrame)
                            F=src.read();
                            
                            %%% Cut off flyback and write frame
                            tgt=Tiff(tgt_name,'a');
                            tgt.setTag('ImageLength',info(1).Height-1);
                            tgt.setTag('ImageWidth', info(1).Width);
                            tgt.setTag('Photometric', Tiff.Photometric.MinIsBlack);
                            tgt.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
                            tgt.setTag('BitsPerSample', 16); % was 8 in example
                            tgt.write(F(1:info(1).Height-1,:))
                            tgt.close();
                            
                            %%% Store flyback info and save to MWorks folder
                            flybacks(iFrame,:)=F(info(1).Height,:);
                        end
                        src.close();
                    end
                end
            end
            
            if self.save_it==true
                %%% Store in MWorks folder
                save(MW_name,'flybacks')
            end
            
            warning on
            
            self.reconstructed=1;
            self.updated=true;
            self.log.append('Pre_processing',t0)
        end
        
        function reconstruct_data(self,varargin)
            error('function not finished...')
            t0=clock;
            nFiles=length(self.files);
            
            warning off
            for iFile=1:nFiles
                tif_name=fullfile(self.dirs.raw_folder,self.files(iFile).name)
                tgt_name=fullfile(self.dirs.rec_folder,self.files(iFile).name);
                
                savec(tgt_name)
                if exist(tif_name,'file')==2
                    info=imfinfo(tif_name);
                    nFrames=length(info);
                    
                    if self.save_it==true
                        flybacks=zeros(nFrames,info(1).Width);
                        src=Tiff(tif_name,'r');
                        
                        for iFrame=1:nFrames
                            src.setDirectory(iFrame)
                            stack=src.read();
                            %stack=uint8(F);
                            size(stack)
                            options.shift =-2;
                            options.DownsamplingRatio = 1;
                            rec = nlwreconstruct(stack,options.shift,options.DownsamplingRatio);
        
                            imshow(rec,[])
                            die
                        end
                        
                    end
                end
            end
            
            self.reconstructed=1;
            self.updated=true;
            self.log.append('Reconstructions',t0)
        end
        
        function register_data(self,varargin)
            %self
            
        end
        
        
        
        
        
        
        
        
        
        
        
        function save_data(self,varargin)
            if self.save_it==true
                if self.updated==true
                    save_name=fullfile(self.dirs.analysis_folder,'dataset.mat');
                    if ~isdir(self.dirs.analysis_folder)
                        mkdir(self.dirs.analysis_folder)
                    end
                    dataset=self;
                    save(save_name,'dataset')
                    fprintf('Data saved to %s!\n',save_name)
                    
                    % Reset flag
                    self.updated=false;
                else
                    disp('No change detected, not saving...')
                end
            else
                disp('Saving disabled...')
            end
        end
        
        function toggle_save(self,varargin)
            if self.save_it==true
                self.save_it=false;
                disp('Disable saving!')
            else
                self.save_it=true;
                disp('Enable saving!')
            end
        end
        
        %function set.save_it(self,varargin)
        %     disp('warning: don''t edit properties directly!!!')
        %end
    end
end