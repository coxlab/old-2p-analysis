clear all
clc

header_script

files=scandir(data_folder,'tif');
nFiles=length(files);
%%

file_name=fullfile(data_folder,files(1).name);
fprintf('Joining datasets in folder of file %s...\n',file_name)
save_name=fullfile(data_folder,'data_analysis',files(1).name);
save_name=strrep(save_name,'tif','mat');
load(save_name,'session_data')

%%% Make sure filenames are relative to data folder on this machine
session_data.rebase(data_root)


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
        session_nr=session_vector(iSession);
        load_name=fullfile(session_data.folder_info.save_folder,files(session_nr).name);
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
        dataset(iClust)=S.join_data_sessions();
        dataset(iClust).cluster_nr=iClust;
        dataset(iClust).session_vector=session_vector;
        dataset(iClust).nFrames=size(dataset(iClust).STIM,1);
    end    
end



%%% Do something with dataset
%% Plot FOV
dataset(2).plot_FOV()

%% look at temporal profile
dataset(2).plot_traces()

%%
dataset(2).RF_analysis(65)

