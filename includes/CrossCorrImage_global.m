function clust=CrossCorrImage_global(varargin)
% BV20150413: goal is to link pixels together that fire together.
% This way, we could average over soma and process pixels to get better
% signal to noise.
%
% !!!! Still under construction !!!!

input=varargin{1};
frames=input.Object;
F=reshape(shiftdim(frames,2),size(frames,3),[]);

sigma=std(F);
TH=prctile(sigma(:),99);
sel=sigma>TH;
MIP=F(:,sel==1);

% Calc distances between them using corr or other measure
CC=corr(MIP);

% Cluster pixels using some clustering method
Z = linkage(CC,'average');
T = cluster(Z,'maxclust',20);

% Construct bwlabel plane
clust=zeros(size(frames,1),size(frames,2));
clust(sel)=T;


