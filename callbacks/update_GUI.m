function update_GUI(varargin)
H=varargin{1};
handles=guidata(H);

%%% Update MIP
MIP=calc_gamma(handles.MIP_raw,handles.MIP_gamma_val);
MIP(MIP(:)<0)=0;
%MIP=imresize(MIP,size(MIP).*[2 1]);
set(handles.subplots(1).h(1),'Cdata',MIP);

set(handles.MIP_selector,'value',handles.MIP_type)

%%% Create list names
if handles.nROI==0
    set(handles.global_properties_table,'String','empty','value',1)
else
    %%% Populate listbox with ROI names
    ROI_nr_vector=cat(2,handles.ROI.ROI_nr);
    handles.ROI_names=strcat({'ROI #'},num2str((ROI_nr_vector).')).';
    set(handles.global_properties_table,'String',handles.ROI_names)
    set(handles.global_properties_table,'value',handles.ROI_selector)
    
    %%% Show properties
    ROI=handles.ROI(handles.ROI_selector);
    
    if ~isempty(ROI.ROI_nr)
        handles.ROI_properties.ROI_nr=ROI.ROI_nr;
        handles.ROI_properties.X_center=ROI.center_coords(1);
        handles.ROI_properties.Y_center=ROI.center_coords(2);
        if isfield(handles,'ellipse_properties')
            handles.ROI_properties.Long_axis=ROI.ellipse_properties.long_axis;
            handles.ROI_properties.Short_axis=ROI.ellipse_properties.short_axis;
            handles.ROI_properties.Radius=mean([ROI.ellipse_properties.long_axis ROI.ellipse_properties.short_axis]);
        end
        
        %%% Show detail
        if 1
            %%            
            %T=handles.MIP(ROI.ROI_rect(2):ROI.ROI_rect(4),ROI.ROI_rect(1):ROI.ROI_rect(3));            
            T=crop_selection(handles.MIP,fliplr(round(ROI.center_coords)),handles.window_size);
            switch 2
                case 0
                    % do nothing
                    disp('Showing unprocessed detail')
                case 1
                    G_size=[3 3];
                    G=bellCurve2(1,G_size/2+1,[1 1],G_size,0);
                    T=convn(T,G);
                    T=T(G_size(1)-1:end-1,G_size(2)-1:end-1);
                case 2
                    %FFT_mask=drawCircle2(size(T,1)/4,[0 0],size(T));                    
                    filters=[.02 .3];
                    FFT_mask_low=drawCircle3(size(T,1)*filters(1),size(T)/2,size(T));
                    FFT_mask_high=drawCircle3(size(T,1)*filters(2),size(T)/2,size(T));
                    FFT_mask=FFT_mask_high-FFT_mask_low;
                    DC=mean(T(:));
                    F=fft2(T-DC);
                    T=real(ifft2(fftshift(fftshift(F).*FFT_mask)))+DC;
            end
            
            fraction=.5;
            T=T/max(T(:));
            
            G=bellCurve2(1,size(T)/2+1,size(T)/3,size(T),0);
            T=T.*G;
            
            N=double(ROI.mask_neuropil)/max(ROI.mask_neuropil(:))*(1-fraction);
            S=double(ROI.mask_soma)/max(ROI.mask_soma(:))*(1-fraction);
            RGB=cat(3,N,T,S);
            RGB(RGB<0)=0;
            RGB=calc_gamma(RGB,1.7);
            RGB=RGB/max(RGB(:));
            handles.detail=RGB;            
            
            set(handles.subplots(2).h(1),'cData',RGB);
            set(handles.subplots(2).fig,'cLim',[min(RGB(:)) max(RGB(:))]);
            %set(handles.subplots(2).fig,'cLim',prctile(RGB(:),[0 99.999]));
        end
        
        %%% Show all ROIs
        %A=cat(1,handles.ROI.coords_MIP_plot);
        %set(handles.subplots(1).p(1),'xData',A(:,1),'yData',A(:,2));        
        for iROI=1:handles.nROI
            A=handles.ROI(iROI).coords_MIP;            
            set(handles.subplots(1).p(10+iROI),'xData',A(:,1),'yData',A(:,2));
        end        
        set(handles.subplots(1).p(2),'xData',ROI.coords_MIP(:,1),'yData',ROI.coords_MIP(:,2));
    else
        handles.nROI=0;
    end
    
    guidata(H,handles)
    writeTable(handles.ROI_properties_table,handles.ROI_properties_table_name)
end
