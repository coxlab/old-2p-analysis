function varargout=crop_selection(varargin)
%function crop_selection(im,coord,patch_rect)

if nargin==0
    %%
    A=rand(100,100);
    patch_rect=[0 0 40 40];    
    rect=OffsetRect(patch_rect,100,50);
    [B,rect_cropped]=crop_selection(A,rect);
    
    imshow(B,[])
    rect_cropped
else
    
    
    if nargin>=1
        im=varargin{1};
    else
        error('Function needs to be called with 2 arguments, 0 given...')
    end
    
    if nargin>=2
        coords=varargin{2};
    else
        error('Function needs to be called with 2 arguments, 1 given...')
    end
    
    
    if nargin>=3
        patch_rect=varargin{3};
        if numel(patch_rect)==2
            patch_rect=[0 0 patch_rect];
        end
    else
        patch_rect=[0 0 40 40];
    end    
    
    %im_rect=[0 0 size(im)];
    src_rect=CenterRectOnPoint(patch_rect,coords(1),coords(2));
    border_size=SizeOfRect(src_rect)/2; % get size of patch
    
    % increase size of source image, to allow selection outside
    im_with_border=centerOnRect(im,size(im)+border_size*2,0);    
    
    % increase coords to keep center in the center
    offset=repmat(border_size,1,2);    
    src_rect_border=src_rect+offset(:)';
    
    % cut out adjusted coords from expanded image
    selected_im=im_with_border(src_rect_border(1)+1:src_rect_border(3),src_rect_border(2)+1:src_rect_border(4));
    
    if nargout>=1
        varargout{1}=selected_im;
    end
    
    if nargout>=2
        varargout{2}=src_rect_border;
    end
    
end