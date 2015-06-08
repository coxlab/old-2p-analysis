function moveFcn_stack(varargin)
figure_handle=varargin{1};
indicator=varargin{3};
target=varargin{4};
frames=varargin{5};

P=get(indicator,'parent');
set(P,'Units','Pixels')

rect=get(P,'position');
current_point=get(figure_handle,'CurrentPoint');
X=current_point(1);
Y=current_point(2);
if between(X,[rect(1) sum(rect([1 3]))])&&between(Y,[rect(2) sum(rect([2 4]))])
    currentPoint=get(P,'CurrentPoint');
    timepoint=floor(currentPoint(1,1));
    
    if between(timepoint,[1 size(frames,3)])
        im=frames(:,:,timepoint);
        set(target,'Cdata',im)
        set(indicator,'Xdata',[timepoint timepoint])
        drawnow
    end
end



