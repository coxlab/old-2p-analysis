function delROI(varargin)
H=varargin{1};
handles=guidata(H);

if handles.nROI>0
    %%% Remove selected ROI
    handles.ROI(handles.ROI_selector)=[];
    
    %%% Remove it from the MIP panel    
    set(handles.subplots(1).p(10+handles.ROI_selector),'xData',[],'yData',[]);
    
    %%% select previous value
    if handles.ROI_selector==1
    else
        handles.ROI_selector=handles.ROI_selector-1;
    end
    set(handles.global_properties_table,'value',handles.ROI_selector)
    
    handles.nROI=length(handles.ROI);
    
    %%% If we just deleted the last one, remove active ROI indicator
    if handles.nROI==0
        set(handles.subplots(1).p(1),'xData',[],'yData',[])
        set(handles.subplots(1).p(2),'xData',[],'yData',[])
        set(handles.subplots(2).h(1),'cData',handles.subplots(2).blank_im);
        
        %%% Clear properties table
        handles.ROI_properties=struct('ROI_nr',[],'X_center',[],'Y_center',[],'Long_axis',[],'Short_axis',[],'Radius',[]);
        guidata(H,handles)
        writeTable(handles.ROI_properties_table,handles.ROI_properties_table_name)
    end
    
    guidata(H,handles)
    
    update_GUI(H)
    
end
