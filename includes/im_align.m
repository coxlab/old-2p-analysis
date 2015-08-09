function varargout=im_align(varargin)
%FOV_shifted=offsetIm(FOV,40,67);
%[CC_max,offset]=im_align(FOV_shifted,FOV)

A=varargin{1};
B=varargin{2};

if nargin>=3&&~isempty(varargin{3})
    use_GPU=varargin{3};
else
    use_GPU=0;
end

if use_GPU==1&&gpuDeviceCount>0
    %%%
    disp('Using GPU acceleration')
    
    CC=gather(normxcorr2(gpuArray(A),gpuArray(B)));
    % get shift coordinates relative to biggest image
    CC_max=max(CC(:));
    [i,j]=find(CC==CC_max);
else
    CC=normxcorr2(A,B);
    if strcmpi(class(CC),'gpuArray')==1
        CC=gather(CC);
    end
    % get shift coordinates relative to biggest image
    CC_max=max(CC(:));
    [i,j]=find(CC==CC_max);
end


peakX=j-size(B,2)/2;
peakY=i-size(B,1)/2;

% shift relative to first image
offset=[peakX-size(B,2)/2 peakY-size(B,1)/2];

varargout{1}=CC_max;
varargout{2}=offset;
