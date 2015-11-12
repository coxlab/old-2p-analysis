function varargout=combine_sparse_masks(varargin)
%function im_2D=combine_sparse_masks(sparse_matrix,field_name)

if nargin>=1
   sparse_matrix=varargin{1};
end

if nargin>=2
   field_name=varargin{2};
end

N=length(sparse_matrix);

M=[];
for i=1:N
    M(:,:,i)=full(sparse_matrix(i).(field_name));
end

im_2d=sum(M,3);
im_2d(im_2d>1)=1;
varargout{1}=im_2d;

