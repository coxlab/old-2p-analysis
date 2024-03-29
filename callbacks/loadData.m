function loadData(varargin)
H=varargin{1};
handles=guidata(H);

%%% load here the separate datafile for the selected session, also save to
%%% this file, leave the session overview untouched!
selected_session=handles.selected_session;
session_data=handles.data_sessions(selected_session);
file_name=session_data.name;
%[folder,file_name]=fileparts(session_data.file_name);
set(handles.session_name,'string',file_name)


loadName=fullfile(handles.data_folder,'data_analysis',file_name);
handles.loadName=loadName;
if exist(loadName,'file')
    load(loadName,'session_data')
else
    % if not found, search for same file relative to current data_root
    fprintf('Looking in location relative to current data root: %s',handles.data_root)
    [f,fn,ext]=fileparts(loadName);
    parts1=strsplit(f,filesep);
    parts2=strsplit(handles.data_root,filesep);
    parts=parts1(length(parts2)+1:end-1);
    folder=strjoin(parts,filesep);
    loadName=fullfile(handles.data_root,folder,'data_analysis',[fn ext]);
    if exist(loadName,'file')
        load(loadName,'session_data')
    else
        error('File not found...')
    end
end

if isfield(session_data,'data')
    nFrames=session_data.data(2);
else
    nFrames=session_data.mov_info.nFrames;
end

%session_data.nFrames=nFrames;
if session_data.is_static_FOV()
    if isfield(session_data,'MIP_std')||isprop(session_data,'MIP_std')
        %%% Get MIP projections if present
        disp('Reloading std MIP')
        if isstruct(session_data.MIP_std)
            handles.MIP_type=3;
            handles.MIP=session_data.MIP_std.data;
            handles.MIP=handles.MIP/max(handles.MIP(:))*256*256;
            handles.MIP_raw=handles.MIP;
            handles.MIP_gamma_val=session_data.MIP_std.gamma_val;
            set(handles.MIP_selector,'value',handles.MIP_type)
        else
            handles.MIP=session_data.MIP_std;
        end
        handles.nFrames=nFrames;
    else
        %%% if not, create them and store to file
        session_file_name=session_data.file_name;
        if exist(session_file_name,'file')==2
            % all is well
        else
            % if not found, search for same file relative to current data_root
            [f,fn,ext]=fileparts(session_file_name);
            parts1=strsplit(f,filesep);
            parts2=strsplit(handles.data_root,filesep);
            parts=parts1(length(parts2)+1:end-1);
            folder=strjoin(parts,filesep);
            session_file_name=fullfile(handles.data_root,folder,[fn ext]);
        end
        info=imfinfo(session_file_name);
        Width=info(1).Width;
        Height=info(1).Height;
        nFrames_max=500;
        nFrames=min([nFrames_max nFrames]);
        tic
        frames=zeros(Height,Width,nFrames);
        for iFrame=1:nFrames
            frame=imread(session_file_name,iFrame,'info',info);
            frames(:,:,iFrame)=double(frame);
        end
        fprintf('Loading frames took: %3.1fs\n',toc);
        
        tic
        session_data.MIP_avg=mean(frames,3);
        session_data.MIP_max=max(frames,[],3);
        session_data.MIP_std=std(frames,[],3);
        handles.MIP=session_data.MIP_std;
        fprintf('Calculating STD took: %3.1fs\n',toc);
    end
    
    %%% show mip
    % make sure image is shown at full size
    handles.subplots(1).fig=subplot(121);
    %handles.MIP_gamma_val=1;
    MIP=calc_gamma(handles.MIP_raw,handles.MIP_gamma_val);
    %MIP=imresize(MIP,size(MIP).*[1 1]);
    %handles.subplots(1).h(1)=imshow(MIP,[0 50]);
    handles.subplots(1).h(1)=imagesc(MIP);
    axis off
    pbaspect([session_data.FOV_info.size_um 1])
    %%% For detail window
    pbaspect(handles.subplots(2).ax,[fliplr(session_data.FOV_info.size_um) 1])
    
    hold on
    handles.subplots(1).p(1)=plot(-1,-1,'r-'); % all ROIs
    axis off % hide plot labels
    
    for iROI=1:handles.nROI_preload
        handles.subplots(1).p(10+iROI)=plot(-1,-1,'r-','buttonDownFcn',{@selectROI,iROI}); % individual ROIs
    end
    handles.subplots(1).p(2)=plot(-1,-1,'c-'); % selected ROI
    hold off
    colormap(handles.green)
    set(handles.subplots(1).fig,'units','normalized','position',[.05 .01 .4 .98])
    set(handles.subplots(1).h(1),'buttonDownFcn',@clickFcnMIP)
    
    %%% Handle ROI definition data
    handles.session_data=session_data;
    handles.ROI_empty=struct('ROI_nr',[],'base_coord',[],'nCoords',0,'coords',[],'ellipse_properties',struct,'ellipse_coords',[],'coords_MIP',[],'coords_MIP_plot',[],'center_coords',[],'ellipse_coords_centered',[],'ROI_rect',[],'mask_soma',[],'mask_neuropil',[],'timeseries_soma',[],'time_series_neuropil',[]);
    handles.usePoly=0;
    
    
    %%%% Load existing ROI definitions
    %session_data.ROI_definitions
    switch 1
        case 1
            try
                ROIs=get_ROI_definitions(session_data,handles.ROI_definition_nr);
                if isempty(ROIs)
                    die
                end
                handles.ROI=ROIs;
                handles.ROI_selector=1;
                handles.nROI=length(handles.ROI);
                fprintf('Reloading %d ROIs\n',handles.nROI)
            catch
                %%% if not, set up an empty structure
                disp('Creating new ROI structure')
                handles.nROI=0;
                handles.ROI_selector=1;
                handles.ROI=handles.ROI_empty;
            end
        case 2
            if isfield(session_data,'ROI_definitions') % should exist from previous step
                if isfield(session_data.ROI_definitions,'ROI') % new
                    if handles.ROI_definition_nr<=length(session_data.ROI_definitions)
                        if isfield(session_data.ROI_definitions(handles.ROI_definition_nr).ROI,'ROI_nr')
                            %%% Get ROI definitions if present
                            handles.ROI=session_data.ROI_definitions(handles.ROI_definition_nr).ROI;
                            handles.ROI_selector=1;
                            handles.nROI=length(handles.ROI);
                            fprintf('Reloading %d ROIs\n',handles.nROI)
                        else
                            %%% if not, set up an empty structure
                            disp('Creating new ROI structure (new)')
                            handles.nROI=0;
                            handles.ROI_selector=1;
                            handles.ROI=handles.ROI_empty;
                        end
                    else
                        %%% if not, set up an empty structure
                        disp('Creating new ROI structure (new user)')
                        handles.nROI=0;
                        handles.ROI_selector=1;
                        handles.ROI=handles.ROI_empty;
                    end
                elseif isfield(session_data.ROI_definitions,'ROI_nr') % old, here for backward compatibility
                    handles.ROI=session_data.ROI_definitions;
                    handles.ROI_selector=1;
                    handles.nROI=length(handles.ROI);
                    fprintf('Reloading %d ROIs (Legacy data...)\n',handles.nROI)
                else
                    %%% if not, set up an empty structure
                    disp('Creating new ROI structure')
                    handles.nROI=0;
                    handles.ROI_selector=1;
                    handles.ROI=handles.ROI_empty;
                end
            end
    end
    handles.ROI_temp=handles.ROI_empty;
    set(handles.global_properties_table,'value',handles.ROI_selector)
    
    guidata(H,handles)
    
    set(handles.subplots(1).p(1),'xData',[],'yData',[])
    set(handles.subplots(1).p(2),'xData',[],'yData',[])
    set(handles.subplots(2).p(1),'xData',[],'yData',[])
    set(handles.subplots(2).p(2),'xData',[],'yData',[])
    set(handles.subplots(2).h(1),'cData',handles.subplots(2).blank_im);
    
    update_GUI(H)
end
