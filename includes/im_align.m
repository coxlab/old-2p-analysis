function varargout=im_align(varargin)
%FOV_shifted=offsetIm(FOV,40,67);
%[CC_max,offset]=im_align(FOV_shifted,FOV)

A=varargin{1};
B=varargin{2};

CC=normxcorr2(A,B);

% get shift coordinates relative to biggest image
CC_max=max(CC(:));
[i,j]=find(CC==CC_max);

peakX=j-size(B,2)/2;
peakY=i-size(B,1)/2;

% shift relative to first image
offset=[peakX-size(B,2)/2 peakY-size(B,1)/2];

varargout{1}=CC_max;
varargout{2}=offset;
