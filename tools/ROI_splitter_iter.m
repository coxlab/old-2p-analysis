clear all
clc

A(:,:,1)=drawCircle2(7,[25 25],[100 100]);
A(:,:,2)=drawCircle2(10,[25 10],[100 100]);
im=sum(A,3)>0;
imshow(im,[])

%%
[decision, result]=split_ROI(im);
%imshow(result,[])

%% 
imshow(,[])

%%

% im_temp=im;
% SE=[0 1 0 ; 1 1 1 ; 0 1 0];
% running=1;
% iter=1;
% 
% while running==1
%     im_temp=imerode(im_temp,SE);
%     CC=bwconncomp(im_temp);
%     N=CC.NumObjects;
%     if 0
%         %%
%         imshow(im_temp,[])
%         title(N)
%     end
%     if N==2
%         disp('2 regions => split') 
%         running=0;
%     end
%     if sum(im_temp(:))==0
%         disp('1 region => don''t split')
%         running=0;
%     end
%     
%     if iter>100
%         running=0;
%     end
% end




