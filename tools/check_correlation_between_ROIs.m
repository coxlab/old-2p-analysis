clear all
clc

header_script

files=scandir(fullfile(data_folder,'data_analysis'),'mat');
nFiles=length(files);

nCols=ceil(sqrt(nFiles));
nRows=ceil(nFiles/nCols);

for iFile=1:nFiles
    
    %save_name=fullfile(data_folder,'data_analysis',files(iFile).name)
    %save_name=strrep(save_name,'tif','mat');
    load_name=fullfile(data_folder,'data_analysis',files(iFile).name);
    load(load_name,'session_data')
    
    if ~session_data.is_static_FOV
        fprintf('Skipping %d',iFile)
    else
        R=corr(session_data.Activity_traces.spike_matrix);
        R(eye(size(R))==1)=0;
        X=squareform(R);
        %%
        subplot(nRows,nCols,iFile)
        %imagesc(triu(R))
        %set(gca,'CLim',[-1 1])
        plot(sort(X,'descend'))
        
        axis([0 length(X) -1 1])
        axis square
        set(gca,'ButtonDownFcn',{@switchFcn,get(gca,'position')})
        drawnow
        
        
        %%
        TH=.2;
        [x,y]=find(triu(R)>TH);
        N=length(x);
        
        % plot ROIs involved
        %length(get_ROI_definitions(session_data))
        %size(session_data.ROI_matrix,2)
        A=combine_sparse_masks(session_data.ROI_matrix(x),'mask_soma');
        B=combine_sparse_masks(session_data.ROI_matrix(y),'mask_soma');
        
        imshow(A-B,[])
        title(sprintf('File nr %d, %d pairs@TH=%3.2f',[iFile N TH]))
        
        if 0
            %% assess individual traces
            iROI=1;
            
            RD=get_ROI_definitions(session_data);
            
            A=session_data.Activity_traces.activity_matrix(:,x(iROI));
            B=session_data.Activity_traces.activity_matrix(:,y(iROI));
            real_corr=corr([A B]);
            subplot(211)
            plot(A)
            hold on
            plot(B+max(A))
            hold off
            title(R(x(iROI),y(iROI)))
            xlabel(real_corr(1,2))
            
            A=combine_sparse_masks(session_data.ROI_matrix(x(iROI)),'mask_soma');
            B=combine_sparse_masks(session_data.ROI_matrix(y(iROI)),'mask_soma');
            subplot(212)
            imshow(A-B,[])
            title([RD(x(iROI)).ROI_nr RD(y(iROI)).ROI_nr])
            
        end
        
    end
end