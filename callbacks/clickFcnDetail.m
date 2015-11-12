function clickFcnDetail(varargin)
H=varargin{1};
handles=guidata(H);

current_point=get(handles.subplots(2).fig,'CurrentPoint');
window_size=handles.window_size(1);
coord=(current_point(1,1:2));
x=coord(1);
y=coord(2);

if isfield(handles,'ROI_temp')
    ROI_temp=handles.ROI_temp;
    if handles.status==0
        
    else
        if between(x,[0 window_size])&&between(y,[0 window_size])
            
            %ROI_selector=handles.ROI_selector;
            if handles.modifiers(1)==0
                %%% save
                ROI_temp.nCoords=ROI_temp.nCoords+1;
                ROI_temp.coords(ROI_temp.nCoords,:)=[x y];
            elseif handles.modifiers(1)==1
                if ROI_temp.nCoords>0
                    %%% remove last added
                    ROI_temp.coords(ROI_temp.nCoords,:)=[];
                    ROI_temp.nCoords=ROI_temp.nCoords-1;
                end
            end
            
            %%% plot definition points
            set(handles.subplots(2).p(1),'xData',ROI_temp.coords(:,1),'yData',ROI_temp.coords(:,2))

            %%% Fit ellipse if not in poly mode
            if ROI_temp.nCoords>5 % only show when enough coordinates are defined
                if handles.usePoly==1
                    poly_coords=ROI_temp.coords;
                    avg_coords=repmat(mean(poly_coords),size(poly_coords,1),1);
                    
                    sort_coords=0;
                    if sort_coords==1
                        poly_coords=poly_coords-avg_coords;
                        % sort in pol space
                        [theta,rho]=cart2pol(poly_coords(:,1),poly_coords(:,2));
                        sorted=sortrows([theta rho],-1);
                        [X,Y]=pol2cart(sorted(:,1),sorted(:,2));
                        poly_coords=[X+avg_coords(:,1) Y+avg_coords(:,2)];                        
                    end                      
                    
                    %%% Add first coordinate at the end to close poly
                    poly_coords=cat(1,poly_coords,poly_coords(1,:));
                    
                    ellipse_coords=poly_coords;
                    offset=avg_coords(1,:)-window_size/2;
                else
                    [ellipse_properties,ellipse_coords]=fit_ellipse(ROI_temp.coords(:,1),ROI_temp.coords(:,2),[]);
                    ellipse_coords=ellipse_coords';
                    offset=[ellipse_properties.X0_in ellipse_properties.Y0_in]-window_size/2;
                    ROI_temp.ellipse_properties=ellipse_properties;                    
                end
                                                
                %%% Fit ellipse on coords
                coords_MIP=ellipse_coords+repmat(ROI_temp.base_coord-window_size/2,size(ellipse_coords,1),1);                                
                                                                
                ROI_temp.ellipse_coords=ellipse_coords;
                ROI_temp.ellipse_coords_centered=ellipse_coords-repmat(offset,size(ellipse_coords,1),1);
                
                ROI_temp.center_coords=ROI_temp.base_coord+offset;
                ROI_temp.ROI_rect=round([ROI_temp.center_coords ROI_temp.center_coords])+[-window_size/2+1 -window_size/2+1 window_size/2 window_size/2];
                
                ROI_temp.coords_MIP=coords_MIP;                
                
                %%% show fit on both panels
                set(handles.subplots(2).p(2),'xData',ellipse_coords(1:end-1,1),'yData',ellipse_coords(1:end-1,2)) % detail
                set(handles.subplots(1).p(2),'xData',coords_MIP(:,1),'yData',coords_MIP(:,2),'color','c') % MIP
            else % remove fit from both panels
                set(handles.subplots(2).p(2),'xData',[],'yData',[])
                set(handles.subplots(1).p(2),'xData',[],'yData',[],'color','c')
            end
                                   
            %%% save
            handles.ROI_temp=ROI_temp;
            guidata(H,handles)
        end
    end
end