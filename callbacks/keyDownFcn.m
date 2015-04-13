function keyDownFcn(varargin)
H=varargin{1};
key_event=varargin{2};
handles=guidata(H);

switch key_event.Key
    % modifiers
    case 'shift'
        handles.modifiers(1)=1;
    case 'control'
        handles.modifiers(2)=1;
    case 'alt'
        handles.modifiers(3)=1;
    case '0'  % command
        handles.modifiers(4)=1;
        
        
        % shortcuts
    case 'return'
        switch handles.status
            case 0
                disp('not doing anything')
            case 1
                
                if isfield(handles,'ROI_temp')
                    ROI_temp=handles.ROI_temp;
                    if ~isempty(ROI_temp.ellipse_coords)
                        disp('Adding ROI')
                        
                        if handles.nROI==0
                            new_nr=1;
                        else
                            new_nr=max(cat(1,handles.ROI.ROI_nr))+1;
                        end
                        ROI_temp.ROI_nr=new_nr;
                        handles.nROI=handles.nROI+1;
                        
                        try
                            %%% Cut out centered ROI and apply region selections
                            T=handles.MIP(ROI_temp.ROI_rect(2):ROI_temp.ROI_rect(4),ROI_temp.ROI_rect(1):ROI_temp.ROI_rect(3));
                        catch
                            error('Centered coordinates fall outside of image...')
                        end
                        
                        %%% smooth T
                        switch 2
                            case 1
                                G_size=[3 3];
                                G=bellCurve2(1,G_size/2+1,[1 1],G_size,0);
                                T=convn(T,G);
                                T=T(G_size(1)-1:end-1,G_size(2)-1:end-1);
                            case 2
                                FFT_mask=drawCircle2(size(T,1)/4,[0 0],size(T));
                                DC=mean(T(:));
                                F=fft2(T-DC);
                                T=real(ifft2(fftshift(fftshift(F).*FFT_mask)))+DC;
                        end
                        
                        %%% Generate mask to isolate soma pixels
                        mask_soma=poly2mask(ROI_temp.ellipse_coords_centered(:,1),ROI_temp.ellipse_coords_centered(:,2),size(T,1),size(T,2));
                        
                        %%% Generate mask to isolate neuropil pixels
                        mask_neuropil=imresize(mask_soma,sqrt(2),'nearest');
                        ext=round((size(mask_neuropil)-handles.window_size)/2);
                        mask_neuropil=mask_neuropil(ext:ext+handles.window_size(1)-1,ext:ext+handles.window_size(2)-1);
                        mask_neuropil=mask_neuropil-mask_soma;
                        
                        %%% Save to main ROI structure
                        ROI_temp.mask_soma=mask_soma;
                        ROI_temp.mask_neuropil=mask_neuropil;
                        handles.ROI(handles.nROI)=ROI_temp;
                        handles.ROI_selector=handles.nROI; % make this ROI the selected one
                        
                        %%% Clear data, ready for next ROI
                        handles.ROI_temp=handles.ROI_empty;
                        
                        %%% next section could be moved to the update gui fcn
                        % we also want to run this after selecting a new ROI from the
                        % list
                        
                        handles.ROI_selector=handles.nROI;
                        set(handles.global_properties_table,'value',handles.ROI_selector)
                        
                        guidata(H,handles)
                        
                        %%% Clear gui
                        set(handles.subplots(1).p(2),'xData',[],'yData',[])
                        set(handles.subplots(2).p(1),'xData',[],'yData',[])
                        set(handles.subplots(2).p(2),'xData',[],'yData',[])
                        
                        update_GUI(H)
                    else
                        disp('not enough data')
                    end
                    
                    % signal job is complete
                    handles.status=0;
                end
            case 2
                disp('editing')
                
                %%% update existing ROI properties based on those in
                %%% ROI_temp
                handles.ROI_selector
                ROI_temp=handles.ROI_temp;
                
                T=handles.MIP(ROI_temp.ROI_rect(2):ROI_temp.ROI_rect(4),ROI_temp.ROI_rect(1):ROI_temp.ROI_rect(3));
                
                %%% Generate mask to isolate soma pixels
                mask_soma=poly2mask(ROI_temp.ellipse_coords_centered(:,1),ROI_temp.ellipse_coords_centered(:,2),size(T,1),size(T,2));
                
                %%% Generate mask to isolate neuropil pixels
                mask_neuropil=imresize(mask_soma,sqrt(2),'nearest');
                ext=round((size(mask_neuropil)-handles.window_size)/2);
                mask_neuropil=mask_neuropil(ext:ext+handles.window_size(1)-1,ext:ext+handles.window_size(2)-1);
                mask_neuropil=mask_neuropil-mask_soma;
                
                %%% Save to main ROI structure
                ROI_temp.mask_soma=mask_soma;
                ROI_temp.mask_neuropil=mask_neuropil;
                handles.ROI(handles.ROI_selector)=ROI_temp;
                
                %%% Clear data, ready for next ROI
                handles.ROI_temp=handles.ROI_empty;
                
                % signal job is complete
                handles.status=0;
                
                %%% Clear gui
                set(handles.subplots(1).p(2),'xData',[],'yData',[])
                set(handles.subplots(2).p(1),'xData',[],'yData',[])
                set(handles.subplots(2).p(2),'xData',[],'yData',[])
                
                guidata(H,handles)
                update_GUI(H)
        end
    case 's'
        if handles.modifiers(4)==1
            saveData(H)
        else
            disp('dummy save')
        end
        
    case 'p'        
        handles.usePoly=1-handles.usePoly;
        if handles.usePoly==1
            disp('Poly mode')
        else
            disp('Ellipse mode')
        end
        
    otherwise % for any other key
        session_data=handles.session_data;
        session_data.ROIs=handles.ROI;
        %session_data
end

guidata(H,handles)