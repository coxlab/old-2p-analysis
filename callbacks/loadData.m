function loadData(varargin)
H=varargin{1};
handles=guidata(H);

%%% load here the separate datafile for the selected session, also save to
%%% this file, leave the session overview untouched!
selected_session=handles.selected_session;
session_data=handles.data_sessions(selected_session);
[folder,file_name]=fileparts(session_data.file_name);
set(handles.session_name,'string',file_name)

loadName=fullfile(folder,'data_analysis',[file_name '.mat']);
handles.loadName=loadName;
if exist(loadName,'file')
    load(loadName,'session_data')
end

nFrames=session_data.data(2);
session_data.nFrames=nFrames;

if isfield(session_data,'MIP_std')
    %%% Get MIP projections if present
    disp('Reloading std MIP')
    if isstruct(session_data.MIP_std)
        handles.MIP_type=3;
        handles.MIP=session_data.MIP_std.data;
        handles.MIP_gamma_val=session_data.MIP_std.gamma_val;
    else
        handles.MIP=session_data.MIP_std;
    end
    handles.nFrames=session_data.data(2);
else
    %%% if not, create them and store to file
    session_file_name=session_data.file_name;
    if exist(session_file_name,'file')==2
        % all is well
    else
        % if not found, search for same file relative to current data_root        
        [f,fn,ext]=fileparts(session_file_name);
        parts=strsplit(filesep,f);
        folder=parts{end};
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
handles.MIP_gamma_val=1;
MIP=calc_gamma(handles.MIP,handles.MIP_gamma_val);
handles.subplots(1).h(1)=imshow(MIP,[]);
hold on
handles.subplots(1).p(1)=plot(-1,-1,'r-'); % all ROIs

handles.nROI_preload=200;
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
handles.ROI_empty=struct('ROI_nr',[],'base_coord',[],'nCoords',0,'coords',[],'ellipse_properties',struct,'ellipse_coords',[],'coords_MIP',[],'coords_MIP_plot',[],'center_coords',[],'ellipse_coords_centered',[],'ROI_rect',[],'mask_soma',[],'mask_neuropil',[]);
handles.usePoly=0;

if isfield(session_data,'ROI_definitions') && isfield(session_data.ROI_definitions,'ROI_nr')
    %%% Get ROI definitions if present
    disp('Reloading ROIs')
    session_data
    handles.ROI=session_data.ROI_definitions;
    handles.ROI_selector=1;
    handles.nROI=length(handles.ROI);
else
    %%% if not, set up an empty structure
    disp('Creating new ROI structure')
    handles.nROI=0;
    handles.ROI_selector=1;    
    handles.ROI=handles.ROI_empty;
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

