function im=drawCircle2(varargin)
% function im=drawCircle(varargin)
% radius - center - imSize

if nargin==0
    help drawCircle2
    return
end

if nargin>=1
    radius=varargin{1};
else
    radius=10;
end
if nargin>=2
    center=varargin{2};
else
    center=[1 1]*radius*2;
end

if nargin>=3
    imSize=varargin{3};
else
    imSize=center*2;
end


% if nargin==0
%     help drawCircle2
% elseif nargin==1
%     radius=varargin{1};
%     center=[0 0];
%     imSize=[radius*2+1 radius*2+1];
% elseif nargin==2
%     radius=varargin{1};
%     center=varargin{2};
%     imSize=[radius*2+1 radius*2+1];
% elseif nargin==3
%     radius=varargin{1};
%     center=varargin{2};
%     imSize=varargin{3};
% else
%     disp('Wrong input count')
% end


% create patch with correct radius
white=1;
%[x,y]=meshgrid(-radius:radius, -radius:radius);
w=imSize/2;
offset=center-w;
[x,y]=meshgrid((-w+1:w)-offset(1), (-w+1:w)-offset(2));
im= white * (x.^2 + y.^2 <= (radius)^2);

% create image of right size
%im=centerOnRect(circle,imSize);

% copy patch to center position in image
%im=offsetIm(im,center(1),center(2),0);

