clear all
clc

header_script
%%
file_name='C:\Users\LBP\Dropbox (coxlab)\2p-datasets\AH03\renders\Maps\extra_stuff\BG_AH03.png';
A_full=double(imread(file_name))/255;

file_name='C:\Users\LBP\Dropbox (coxlab)\2p-datasets\AH03\renders\Maps\extra_stuff\Cesar\16bitsurf.tif';
B_full=double(imread(file_name))/255;

%% pre-processing
switch 0
    case 0
        A=A_full(:,:,2);
        B=B_full;
    case 1
        ws=20;
        C=0/100;
        tm=0;
        A=-adaptivethreshold(A_full(:,:,2),ws,C,tm);
        B=-adaptivethreshold(B_full,ws,C,tm);
end


%%
subplot(221)
imshow(A,[])
hold on
p(1)=plot(0,0,'m*');
p(6)=plot(0,0,'b.');
hold off
subplot(222)
imshow(B,[])
hold on
p(2)=plot(0,0,'m*');
hold off
axis equal

subplot(2,2,[3 4])
p(3)=plot(0,0,'bs');
hold on
p(4)=plot(0,0,'ro');
p(5)=plot(0,0,'go');
hold off
t=title('');
axis ij
axis square

colormap(green)
switch 2
    case 1 %% normxcorr
        %%
        [CC_max,offset]=im_align(B,A);
        
        subplot(121)
        imshow(A,[])
        hold on
        plot(size(A,1)/2+offset(2),size(A,2)/2+offset(1),'m*')
        hold off
    case 2 % manual
        collecting=1;
        current_pair=1;
        current_plot=1;
        data=struct;
        while collecting==1
            subplot(2,2,current_plot)
            switch current_plot
                case 1
                    [x,y,button]=ginput(1);
                    if button==1
                        data(current_pair).left_coords=[x y];
                        
                        current_plot=2;
                    end
                case 2
                    [x,y,button]=ginput(1);
                    if button==1
                        data(current_pair).right_coords=[x y];
                        
                        current_plot=1;
                        current_pair=current_pair+1;
                    end
            end
            
            if isfield(data,'left_coords')
                L=cat(1,data.left_coords);
                set(p(1),'Xdata',L(:,1),'Ydata',L(:,2))
                set(p(3),'Xdata',L(:,1),'Ydata',L(:,2))
            end
            if isfield(data,'right_coords')
                R=cat(1,data.right_coords);
                set(p(2),'Xdata',R(:,1),'Ydata',R(:,2))
                set(p(4),'Xdata',R(:,1),'Ydata',R(:,2))
                set(p(6),'Xdata',R(:,1),'Ydata',R(:,2))
            end
            
            if current_pair>=3&&current_plot==1
                %%
                [D, Z, TRANSFORM]=procrustes(L,R);
                set(p(5),'Xdata',Z(:,1),'Ydata',Z(:,2))
                set(t,'string',sprintf('%3.2f',D))
                
                Z_reconstruct = TRANSFORM.b * R * TRANSFORM.T + TRANSFORM.c;
            end
            
            drawnow
            
            if ismember(GetKey,[13 27])
                collecting=0;
            end
        end
end


%%
L=cat(1,data.left_coords);
R=cat(1,data.right_coords);

[D, Z, TRANSFORM]=procrustes(L,R);

%%

switch 2
    case 1 % image way
        B_scale=imresize(B,TRANSFORM.b);
        %B_sized=centerOnRect(B_scale,size(A),128);
        %offset=mean(TRANSFORM.c,1)-size(B)/2;
        offset=[0 0];
        B_offset=offsetIm(B_scale,offset(1),offset(2));
        B_combined=A*.5+B_offset/256;
    case 2 % rect way        
        scale=TRANSFORM.b*.95;
        offset=mean(TRANSFORM.c,1)*1;
        B_rect=[0 0 size(B)];
        B_rect=ScaleRect(B_rect,scale,scale);
        B_rect=OffsetRect(B_rect,offset(2),-offset(1));
        B_rect=(B_rect);
        
        B_scaled=imresize(B,scale)/256;
        B_scaled=B_scaled(1:end-10,:);
        
        B_combined=A;
        patch=B_combined(B_rect(1):B_rect(1)+size(B_scaled,1)-1,B_rect(2):B_rect(2)+size(B_scaled,2)-1);
        size(patch)
        size(B_scaled)
        patch=patch+B_scaled;
        B_combined(B_rect(1):B_rect(1)+size(B_scaled,1)-1,B_rect(2):B_rect(2)+size(B_scaled,2)-1)=patch;
        
end

subplot(121)
imshow(A,[])
hold on
plot(B_rect([2 4]),B_rect([1 3]),'r')
hold off
subplot(122)
imshow(B_combined,[])
colormap(green)


