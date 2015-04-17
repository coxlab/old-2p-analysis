clear all
clc

data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-04-15_AF11';
loadName=fullfile(data_folder,'data_analysis','20150415_AF11_RM_004.mat');

if exist(loadName,'file')
    load(loadName,'session_data')
    
    %%
    stimulus_matrix_ext=session_data.stimulus_matrix_ext;
    activity_matrix=session_data.activity_matrix;
    
    for iROI=1
        condition_vector=stimulus_matrix_ext(:,5);
        conditions=unique(condition_vector(condition_vector>0));
        nConditions=length(conditions);
        for iCond=1%:nConditions
            % select all trials for this condition
            condition_nr=conditions(iCond);
            
            sel=stimulus_matrix_ext(:,5)==condition_nr;
            condition_data=stimulus_matrix_ext(sel,:);
            nFrames=size(activity_matrix,1);
            T=(1:nFrames)/3;
            plot(T,activity_matrix(:,[1:end]))
            hold on
            %plot(condition_data(:,1),condition_data(:,6:7))
            plot(stimulus_matrix_ext(:,1),stimulus_matrix_ext(:,5))
            hold off
            axis([0 500 -50 100])
        end
    end
    
end