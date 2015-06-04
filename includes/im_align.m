function varargout=im_align(varargin)

A=varargin{1};
B=varargin{2};

CC=normxcorr2(B,A);

% get shift coordinates relative to biggest image
CC_max=max(CC(:));
[i,j]=find(CC==CC_max);

peakX=j-size(B,2)/2+1;
peakY=i-size(B,1)/2+1;

% shift relative to first image
offset=[peakX-size(B,2)/2 peakY-size(B,1)/2]-[1 0];

varargout{1}=CC_max;
varargout{2}=offset;
