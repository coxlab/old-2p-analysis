function varargout=split_ROI(varargin)
if nargin>=1
    im_temp=varargin{1};
else
    error('Function needs at least 1 input')
end
if nargin>=2
    SE=varargin{2};
else
    SE=[0 1 0 ; 1 1 1 ; 0 1 0];
    SE=strel('disk',3,0);
end

%%% Iterative erode until 1 or more points remain
running=1;
iter=1;
decision=[];
while running==1
    %%% Erode
    im_prev=im_temp;
    im_temp=imerode(im_temp,SE);
    
    %%% Check number of components
    CC=bwconncomp(im_temp);
    N=CC.NumObjects;
    
    if 1
        %% Show change over iterations
        imshow(im_temp,[])
        title(N)
    end
    
    if sum(im_temp(:))==0
        running=0;
        decision=0;
        disp('1 region => don''t split')
    elseif N==2
        running=0;
        decision=1;
        disp('2 regions => split')
        
        %imshow(im_prev-im_temp,[])
        
        im_label=bwlabel(im_temp);
        for iROI=1:2
            I=im_label==iROI;            
            for i=1:iter+1
                I=imdilate(I,SE);
            end
            result(:,:,iROI)=I;
            %imshow(,[])
        end
        
        %%% Calc overlap region
        im_out=sum(result,3);
        overlap_px=im_out==2;
        
        %%% Remove pixels that overlap with other region
        for iROI=1:2
            result(:,:,iROI)=(result(:,:,iROI)-overlap_px);            
            I(:,:,iROI)=varargin{1}.*result(:,:,iROI)*iROI;
        end
        im_out=sum(I,3);
        
        im_out=kmeans(varargin{1},2)
        imshow(im_out,[]) 
        
        

%         overlap_px=imdilate(overlap_px,SE);
%         overlap_px=imdilate(overlap_px,SE);

%         im_out=sum(result,3);
%         res=varargin{1}-overlap_px;
%         imshow(res,[])
    end
    
    
                
    if iter>100
        running=0;
    end
    iter=iter+1;
end

if nargout>=1
    varargout{1}=decision;
end
if nargout>=2
    varargout{2}=im_out;
end