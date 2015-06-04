function out=centerOnRect(img,varargin)
% This function creates a new image of a certain size or original size +
% border (when given a scalar input). The original image will be centered
% on the new one, and borders will be cropped or padded with zeros
% to obtain requested size image. With only one input, the images will be
% put on a square field

% get original size
imSize=size(img);

if nargin==1
    newSize=[max(imSize) max(imSize)];
elseif nargin>=2
    temp=varargin{1};
    if numel(temp)==1
        newSize=[max(imSize) max(imSize)]+temp;
    elseif numel(temp)==2
        newSize=temp;
    end
end

if nargin<3
    bg=img(1,1);
else
    bg=varargin{2};
end

newSize=ceil(newSize);

sideLength=ceil(max([imSize newSize]));
diffRows=(sideLength-imSize(1))/2;
diffCols=(sideLength-imSize(2))/2;
out_temp=zeros(sideLength)+bg;
out_temp(floor(diffRows)+1:end-ceil(diffRows),floor(diffCols)+1:end-ceil(diffCols))=img;

diffRows=(sideLength-newSize(1))/2;
diffCols=(sideLength-newSize(2))/2;
out=out_temp(floor(diffRows)+1:end-ceil(diffRows),floor(diffCols)+1:end-ceil(diffCols));