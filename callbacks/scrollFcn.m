function scrollFcn(varargin)
H=varargin{1};
key_event=varargin{2};
handles=guidata(H);

if isfield(handles,'detail')
    if handles.status==0
        %%% Scroll through ROIs
        selector_change=key_event.VerticalScrollCount;
        ROI_selector=handles.ROI_selector+selector_change;
        
        if between(ROI_selector,[1 handles.nROI+1])
            handles.ROI_selector=ROI_selector;
            guidata(H,handles)            
            update_GUI(H)
        end
        
    else
        delta_gamma=-key_event.VerticalScrollCount/100;
        handles.detail_gamma_val=handles.detail_gamma_val+delta_gamma;
        if handles.detail_gamma_val>4
            handles.detail_gamma_val=4;
        end
        if handles.detail_gamma_val<1e-6
            handles.detail_gamma_val=1e-6;
        end
        
        D=calc_gamma(handles.detail,handles.detail_gamma_val);
        set(handles.subplots(2).h(1),'cData',D);
        set(handles.subplots(2).fig,'cLim',[min(D(:)) max(D(:))]);
    end
    
end
guidata(H,handles)

