function image=offsetIm(im,X,Y,varargin)

X=round(X);
Y=round(Y);

if nargin<4
    bg=mean(im(:));
else
    bg=varargin{1};
end
    
image=zeros(size(im))+bg;
if X>=0&&Y>=0
    srcRect=RectOfMatrix(im)+[0 0 -abs(X) -abs(Y)];
elseif X>=0&&Y<0
    srcRect=RectOfMatrix(im)+[0 abs(Y) -abs(X) 0];
elseif X<0&&Y>=0
    srcRect=RectOfMatrix(im)+[abs(X) 0 0 -abs(Y)];
elseif X<0&&Y<0
    srcRect=RectOfMatrix(im)+[abs(X) abs(Y) 0 0];
end
dstRect=round(OffsetRect(srcRect,X,Y));
image(dstRect(2)+1:dstRect(4),dstRect(1)+1:dstRect(3))=im(srcRect(2)+1:srcRect(4),srcRect(1)+1:srcRect(3));
