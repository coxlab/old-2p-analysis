%% Join datasets

% check folder for joinable sessions: close in space (and time)
[clusters,cluster_vector,nClusters,files]=session_data.find_FOV_clusters(session_data.folder_info.save_folder);

% if close in space, but not in time, link files but call it a different
% timepoint, so you can look at evolution and/or stability

for iClust=1:nClusters
    cluster_nr=cluster_vector(iClust);
    session_vector=find(clusters==cluster_nr);
    nSessions=length(session_vector);
    
    %% Load all session into one compound object
    clear S
    for iSession=1:nSessions
        load_name=fullfile(session_data.folder_info.save_folder,files(iSession).name);
        load(load_name,'session_data')
        S(iSession)=session_data;
    end
    
    %% Check whether same number of ROIs are defined in all sessions
    ROI_counts=S.get_ROI_counts(ROI_definition_nr);
    if all(ROI_counts==ROI_counts(1))
        nROIs=mean(ROI_counts);
    else
        error('Number of ROIs has to be identical for all sessions');
    end
    
    
    %% Verify distances are closely matched over sessions
    distances=S.get_ROI_distances(ROI_definition_nr);
    if any(mean(distances)>10)
        error('Large shift in ROI centers detected...')
    else
        % all good: take the ROI_definitions of one of the sessions
        %ROI_definitions=S(1).ROI_definitions(ROI_definition_nr).ROI;
        
        % Then join them 
        % leave out last trial
        dataset=S.join_data_sessions(ROI_definition_nr);
        dataset.cluster_nr=iClust;
        dataset.session_vector=session_vector;
        dataset.nFrames=size(dataset.STIM,1);
    end
        
    
    %% Do something with dataset
    dataset.plot_traces()
end