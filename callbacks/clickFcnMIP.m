function clickFcnMIP(varargin)
H=varargin{1};
handles=guidata(H);

if isfield(handles,'session_data')
    
    session_data=handles.session_data;
    
    %%% Get general movie properties
    if isfield(session_data,'data')
        Width=session_data.data(3);
        Height=session_data.data(4);
    else
        Width=session_data.mov_info.Width;
        Height=session_data.mov_info.Height;
    end
    
    %%% Read point user click on
    current_point=get(handles.subplots(1).fig,'CurrentPoint');
    coord=round(current_point(1,1:2));
    
    % if on ROI, select this ROI
    nROI=handles.nROI;
    in_polygon_vector=zeros(nROI,1);
    for iROI=1:nROI
        if inpolygon(coord(1),coord(2),handles.ROI(iROI).coords_MIP(:,1),handles.ROI(iROI).coords_MIP(:,2))
            in_polygon_vector(iROI)=1;
        end
    end
    
    if any(in_polygon_vector)
        % select ROI
        selectROI(H,find(in_polygon_vector==1))
    else
        % otherwise, start creation of new ROI
        handles.status=1;
        ext=handles.window_size(1)/2;
        x=coord(2);
        y=coord(1);
        if between(x,[ext Height-ext])&&between(y,[ext Width-ext])
            %%% Select corresponding part of MIP
            detail=handles.MIP(x-ext+1:x+ext,y-ext+1:y+ext);
            handles.detail_gamma_val=1;
            
            %%% Apply slight blur to get rid of pixelation
            switch 2
                case 0
                    disp('Showing unprocessed detail')
                case 1
                    G_size=[3 3];
                    G=bellCurve2(1,G_size/2+1,[1 1],G_size,0);
                    detail=convn(detail,G);
                    detail=detail(G_size(1)-1:end-1,G_size(2)-1:end-1);
                case 2
                    FFT_mask=drawCircle2(size(detail,1)/4,[0 0],size(detail));
                    DC=mean(detail(:));
                    F=fft2(detail-DC);
                    detail=real(ifft2(fftshift(fftshift(F).*FFT_mask)))+DC;
            end
            switch 2
                case 1
                    %%% Save
                    handles.ROI_selector=handles.nROI+1;
                    handles.ROI(handles.ROI_selector).base_coord=coord;
                    handles.ROI(handles.ROI_selector).nCoords=0;
                    handles.ROI(handles.ROI_selector).coords=[];
                    handles.detail=detail;
                case 2
                    handles.detail=detail;
                    ROI_temp=handles.ROI_empty;
                    ROI_temp.base_coord=coord;
                    handles.ROI_temp=ROI_temp;
            end
            
            %%% Get activity trace for whole window
            %handles.data_sessions
            
            %%% Show
            set(handles.subplots(1).p(2),'xData',[],'yData',[])
            set(handles.subplots(2).fig,'cLim',[min(handles.detail(:)) max(handles.detail(:))]);
            set(handles.subplots(2).h(1),'cData',handles.detail);
            set(handles.subplots(2).p(1),'xData',[],'yData',[])
            set(handles.subplots(2).p(2),'xData',[],'yData',[])
        else
            
        end
        guidata(H,handles)
    end    
end