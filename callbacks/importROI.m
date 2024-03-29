function importROI(varargin)
H=varargin{1};
if nargin>=3&&~isempty(varargin{3})
    source=varargin{3};
else
    source=1;
end

handles=guidata(H);

switch source
    case 1 % another session
        %%% Get file to load ROIs from
        cd(fullfile(handles.data_folder,'data_analysis'))
        [filename,pathname]=uigetfile('.mat');
        loadName=fullfile(pathname,filename);
        
        %%% Load data from file
        load(loadName,'session_data')
        
        %%% Load ROI definitions from selected session_data
        ROI=get_ROI_definitions(session_data,handles.ROI_definition_nr);
        nROI=length(ROI);
        
        %%% implement shift to coords, based on offset between MIP images
        im1=handles.MIP; % MIP from current file
        if isfield(session_data,'MIP_std')
            if isfield(session_data.MIP_std,'data')
                im2=session_data.MIP_std.data; % temp MIP from other, existing file
            else
                im2=session_data.MIP_std; % for older version old
            end
        elseif isprop(session_data,'MIP_std')
            im2=session_data.MIP_std.data; % temp MIP from other, existing file
        else
            session_data
            die
        end
        
        %%% Check GUI checkbox whether to align automatically or not
        if get(handles.auto_align,'value')==1
            [CC_max,offset]=im_align(im1,im2);
            
            if CC_max>.90
            else
                disp('Max cross-correlation was too low, no alignment was applied. Are you sure these fields are the same?')
                offset=[0 0];
            end
%             if numel(im1)>numel(im2)
%                 A=im1;
%                 template=im2;
%                 ref=1;
%             else
%                 A=im2;
%                 template=im1;
%                 ref=2;
%             end
%             CC=normxcorr2(template,A);
%             
%             % get shift coordinates relative to biggest image
%             CC_max=max(CC(:));
%             [i,j]=find(CC==CC_max);
%             
%             peakX=j-size(template,2)/2+1;
%             peakY=i-size(template,1)/2+1;
%             
%             if ref==1 % make shift relative to first image
%                 offset=-[peakX-size(template,2)/2 peakY-size(template,1)/2]-[1 0];
%             else
%                 offset=[peakX-size(template,2)/2 peakY-size(template,1)/2]-[1 0];
%             end
        else
            offset=[0 0];
        end
        
    case 2 % auto ROIs
        ROI=get_ROI_definitions(handles.session_data,1);
        nROI=length(ROI);
        
        offset=[0 0];
        loadName=handles.session_data.save_name;
end

% move all coordinates with the offset found
if nROI==0
    disp('No ROIs found')
else
    for iROI=1:nROI
        ROI(iROI).base_coord=ROI(iROI).base_coord-offset;
        ROI(iROI).coords=ROI(iROI).coords-repmat(offset,size(ROI(iROI).coords,1),1);
        ROI(iROI).ellipse_coords=ROI(iROI).ellipse_coords-repmat(offset,size(ROI(iROI).ellipse_coords,1),1);
        ROI(iROI).coords_MIP=ROI(iROI).coords_MIP-repmat(offset,size(ROI(iROI).coords_MIP,1),1);
        %ROI(iROI).coords_MIP_plot=ROI(iROI).coords_MIP_plot-repmat(offset,size(ROI(iROI).coords_MIP_plot,1),1);
        ROI(iROI).center_coords=ROI(iROI).center_coords-repmat(offset,size(ROI(iROI).center_coords,1),1);
        ROI(iROI).ellipse_coords_centered=ROI(iROI).ellipse_coords_centered-repmat(offset,size(ROI(iROI).ellipse_coords_centered,1),1);
        ROI(iROI).ROI_rect=ROI(iROI).ROI_rect-repmat(offset,size(ROI(iROI).ROI_rect,1),2);
    end
    fprintf('Shifted all coordinates by x=%d and y=%d \n',offset)
    
    %%% Check if coords are still valid after the shift
    if isfield(handles.session_data,'data')
        Height=handles.session_data.data(4);
        Width=handles.session_data.data(3);
    else
        Height=handles.session_data.mov_info.Height;
        Width=handles.session_data.mov_info.Width;
    end
    
    if 0
        remove_vector=zeros(nROI,1);
        for iROI=1:nROI
            ROI_rect=ROI(iROI).ROI_rect;
            %if any(~between(ROI_rect([1 3]),[1 Height]))&&any(~between(ROI_rect([2 4]),[1 Width]))
            if any(~between(ROI_rect([1 3]),[1 Width]))||any(~between(ROI_rect([2 4]),[1 Height]))
                %disp('ROI is invalid after shift, removing')
                remove_vector(iROI)=1;
            else
                %disp('ROI still valid')
            end
        end
        ROI(remove_vector==1)=[];
        nRemoved=sum(remove_vector);
        if nRemoved>0
            fprintf('Removed %d ROI that became invalid after shift \n',nRemoved)
        end
    end
    
    %%% Import ROIs from other file
    handles.ROI=ROI;
    handles.ROI_selector=1;
    handles.nROI=length(handles.ROI);
    handles.ROI_temp=handles.ROI_empty;
    
    set(handles.global_properties_table,'value',handles.ROI_selector)
    fprintf('Imported %d out of %d ROIs from %s \n',[handles.nROI nROI loadName])
    
    guidata(H,handles)
    
    update_GUI(H)
end
