function step03_ROI_GUI(varargin)
% BV20150309: This function will take a series of Z-projections from a movie and allow
% the user with a bunch of cool tools to select all ROIs contained within
% the field.
% Press 'n' to start new ROI definition, first click points to the approx
% center and all subsequent clicks are used to the fit the outline.
% Ellipse or poly in case of weird shapes.
% Links between files can be made and files can be chopped up and linked in
% case of shifts during the session.

% Layout:
% - main window to select ROIs
% - detail window to finetune selected ROI
%        contrast can be changed, independent from the image
% - list view so you can select and edit ROIs: need to get back into tables
% - list view of ROI properties/details

% 2DO:
% V add option to load ROIs from other file
%   It is recommended to load existing ROI definitions for overlapping FOVs
%   as to maintain the cell identity over sessions/days. That is why step02
%   is about visualizing the FOVs and finding overlapping sessions
% V Create layer on top of images to enable clicking on the image even when
%   lines are plotted on top of it. Did something similar for cluster cutting
%   GUI, need to make sure both axes are identical size at all times.
% V add gamma correct for main MIP window, to bring out more detail

% BV20150408: compare and use this as ground truth for developing automated
% methods

% BV20150409: save gamma val for MIP if it is set and reload
% BV20150413: allow switching between MIPs


%%% Add subfolders
path_dir=fileparts(mfilename('fullpath'));
addpath(genpath(path_dir))

%%% Set data_root
header_script

handles.data_root=data_root;
handles.green=green;
handles.window_size=[40 40];

%%% Set up GUI
handles.figure1=figure(999);
clf
set(handles.figure1,'units','normalized','position',[0.60 0.25 0.35 0.5],'resize','on','menubar','none','NumberTitle','Off','Name','ROI_GUI: ROI extraction from 2p data')
mainPanel=uipanel(handles.figure1,'units','normalized','position',[.01 .01 .98 .98]);

graphPanel=uipanel(mainPanel,'units','normalized','position',[.01 .52 .98 .47]);
propPanel=uipanel(mainPanel,'units','normalized','position',[.01 .10 .98 .40]);
controlPanel=uipanel(mainPanel,'units','normalized','position',[.01 .01 .98 .08]);

%%% Enable global keyboard modifiers/shortcuts
set(handles.figure1,'WindowKeyPressFcn',@keyDownFcn)
set(handles.figure1,'WindowKeyReleaseFcn',@keyUpFcn)
set(handles.figure1,'WindowScrollWheelFcn',@scrollFcn)

handles.modifiers=zeros(3,1);
handles.MIP_type=0;
handles.MIP_gamma_val=1;
handles.detail_gamma_val=1;
handles.status=0;

%%% Prepare image placeholder for MIP
blank_im=zeros(512,512);
handles.subplots(1).fig=subplot(1,2,1,'Parent',graphPanel);
handles.subplots(1).h(1)=imshow(blank_im);
% add gamma value buttons
uicontrol(graphPanel,'Style','pushbutton','units','normalized','position',[.01 .35 .04 .1],'string','^','callback',{@MIP_gamma,'up'})
uicontrol(graphPanel,'Style','pushbutton','units','normalized','position',[.01 .25 .04 .1],'string','v','callback',{@MIP_gamma,'down'})
uicontrol(graphPanel,'Style','popupmenu','units','normalized','position',[.01 .1 .08 .01],'string',{'AVG','MAX','STD','CC'},'callback',@select_MIP)

%%% Prepare image placeholder for detail blow-up
blank_im=zeros(handles.window_size);
handles.subplots(2).fig=subplot(1,2,2,'Parent',graphPanel);
handles.subplots(2).h(1)=imshow(blank_im);
handles.subplots(2).blank_im=blank_im;
hold on
handles.subplots(2).p(1)=plot(-1,-1,'r*'); % marker
handles.subplots(2).p(2)=plot(-1,-1,'c-','buttonDownFcn',@clickFcnDetail); % line fit
hold off
colormap(handles.green)
set(handles.subplots(2).fig,'units','normalized','position',[.55 .01 .4 .98])
set(handles.subplots(2).h(1),'buttonDownFcn',@clickFcnDetail)

%%% Prepare data structures
%global_properties=struct('ROI_nr',0,'Shape','Poly');
ROI_properties=struct('ROI_nr',[],'X_center',[],'Y_center',[],'Long_axis',[],'Short_axis',[],'Radius',[]);

%handles.global_properties=global_properties;
handles.ROI_properties=ROI_properties;
handles.selected_session=0;
handles.ROI_selector=0;

%%% Prepare selection list box
handles.global_properties_table=uicontrol(propPanel,'Style','listbox','Units','Normalized','Position',[.01 .01 .25 .98]);
set(handles.global_properties_table,'String','','Callback',@listCallback)

%%% Prepare properties table
handles.ROI_properties_table=uitable(propPanel,'Tag','input_table','Units','Normalized','Position',[.51 .01 .48 .98]);
handles.ROI_properties_table_name='ROI_properties';
set(handles.ROI_properties_table,'ColumnName',{'Parameter' 'Value'},'ColumnWidth',{150 200},'ColumnEditable',[false true],'ColumnFormat',{'char' ''},'RowName',[],'CellEditCallback',{@readTable,handles.ROI_properties_table_name})

%%% Import ROIs from other file
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.30 .9 .1 .1],'string','Import ROIs','callback',@importROI)

% All ROI move buttons
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.325 .8 .05 .1],'string','^','callback',{@shiftROI,'all','up'})
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.295 .7 .05 .1],'string','<','callback',{@shiftROI,'all','left'})
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.355 .7 .05 .1],'string','>','callback',{@shiftROI,'all','right'})
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.325 .6 .05 .1],'string','v','callback',{@shiftROI,'all','down'})

% Edit delete
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.30 .5 .1 .1],'string','Edit ROI','callback',@editROI)
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.30 .4 .1 .1],'string','Delete ROI','callback',@delROI)

% Current ROI move buttons
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.325 .3 .05 .1],'string','^','callback',{@shiftROI,'single','up'})
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.295 .2 .05 .1],'string','<','callback',{@shiftROI,'single','left'})
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.355 .2 .05 .1],'string','>','callback',{@shiftROI,'single','right'})
uicontrol(propPanel,'Style','pushbutton','units','normalized','position',[.325 .1 .05 .1],'string','v','callback',{@shiftROI,'single','down'})



%%% Create buttons
hor_spacing=.10;
uicontrol(controlPanel,'Style','pushbutton','units','normalized','position',[.02 .25 .1 .5],'string','Set dir','callback',@setDir)
handles.session_selector_ui=uicontrol(controlPanel,'Style','popupmenu','units','normalized','position',[.02+hor_spacing*1 .2 .15 .5],'string','0','callback',@selectSession);
uicontrol(controlPanel,'Style','pushbutton','units','normalized','position',[.02+hor_spacing*2+.05 .25 .1 .5],'string','Reload Data','callback',@loadData)
handles.session_name=uicontrol(controlPanel,'Style','edit','units','normalized','position',[.02+hor_spacing*3+.06 .26 .15 .5],'string','No session');
uicontrol(controlPanel,'Style','pushbutton','units','normalized','position',[.02+hor_spacing*5+.05 .25 .1 .5],'string','Save Data','callback',@saveData)


%%% Store current handles
guidata(handles.figure1,handles)

%%% Fill tables with handles
writeTable(handles.ROI_properties_table,handles.ROI_properties_table_name)
