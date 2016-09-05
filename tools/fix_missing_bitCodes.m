% depends on
% - get_MWorks_info(), which needs an (converted) *.mwk file
% - the existance of frame_info.bitCode_vector, extracted directly from the flyback line

clear all
clc

nBits=4;
are_you_sure=true;
batch_files=1;
scenario=1;
switch scenario
    case 1
        data_root='/Users/benvermaercke/Dropbox (coxlab)/2p-data/';
    case 2
        data_root='/Users/benvermaercke/Desktop/coxlab_stuff/data_analysis';
end


switch batch_files
    case 1
        f=uigetdir(data_root,'Select your datafile');
        switch scenario
            case 1
                data_folder=fullfile(f,'data_analysis');
            case 2
                data_folder=f;
        end
        files=scandir(data_folder,'mat');
        nFiles=length(files);
    case 0
        [fn,data_folder]=uigetfile('*.mat','Select your datafile',data_root);
        files(1).name=fn;
        nFiles=1;
end

%%

for iFile=1:nFiles
    load_name=fullfile(data_folder,files(iFile).name);
    
    try
        if isempty(strfind(load_name,'dataset'))&&isempty(strfind(load_name,'overview'))
            load(load_name)
            
            if session_data.is_static_FOV==1
                fprintf('>> Checking file: %s \n',load_name)
                
                %session_data.folder_info.data_folder
                session_data.rebase(data_root)
                %session_data.folder_info.data_folder
                
                scim_bitCodes=mode(cat(2,session_data.frame_info.bitCode_vector))';
                bit_codes_unique=unique(scim_bitCodes);
                
                if length(bit_codes_unique)==2^nBits
                    disp('Good news: no sign of inactive bitcodes!')
                elseif length(bit_codes_unique)==1
                    error('X Zero active bits found, ignoring file...')
                else
                    % find less frequent bitCodes
                    possible_bitCodes=0:2^nBits-1;
                    freqs=hist(scim_bitCodes,possible_bitCodes);
                    sel=(freqs>max(freqs)/10);
                    
                    %%
                    pattern_matrix=dec2bin(possible_bitCodes(sel));
                    
                    inactive_bits=zeros(nBits,1);
                    for bit=1:nBits
                        bit_values=pattern_matrix(:,end-bit+1);
                        if std(bit_values)==0
                            fprintf('>>> Bit %d was inactive... \n',bit)
                            inactive_bits(end-bit+1)=1;
                        end
                    end
                    nInactive_bits=sum(inactive_bits==1);
                    
                    if nInactive_bits==nBits
                        error('X Zero active bits found, ignoring file...')
                    end
                    
                    if nInactive_bits==0
                        disp('! No bit was found to be inactive all the time...')
                    end
                end
                
                %% Load all stim update events
                session_data.bitCodes.MWorks_info=[];
                try
                    session_data.get_MWorks_info()
                catch
                    %session_data.folder_info.data_folder=data_root;
                    %session_data.get_MWorks_info('20160321_JR005W_rsvp.mwk')
                end
                MWorks_info=session_data.bitCodes.MWorks_info;
                T=double(cat(1,MWorks_info.timestamp));
                MWorks_bitCodes=[cat(1,MWorks_info.event_id) T cat(1,MWorks_info.bit_code) [diff(T)/1e6;3]];
                
                
                %% resample MWorks bitCodes to scim frame_rate,
                frame_rate=session_data.mov_info.frame_rate;
                nFrames_per_bitCode=round(MWorks_bitCodes(:,4)*frame_rate);
                MWork_bitCodes_interpolated=[0 0];
                for iBC=1:size(MWorks_bitCodes,1)
                    MWork_bitCodes_interpolated=cat(1,MWork_bitCodes_interpolated,repmat(MWorks_bitCodes(iBC,[1 3]),nFrames_per_bitCode(iBC),1));
                end
                
                %% use this for cc calculation, will tolerate bit-loss
                cc=normxcorr2(scim_bitCodes,MWork_bitCodes_interpolated(:,2));
                [max_corr,loc]=max(cc);
                if max_corr<1
                    max_corr
                end
                offset=loc-length(scim_bitCodes);
                new_vector=MWork_bitCodes_interpolated(offset+1:offset+length(scim_bitCodes),:);
                
                if 0
                    %%
                    plot([new_vector(:,2) scim_bitCodes])
                end
                
                %% now from this new_vector, we can recall all stim properties
                condition_vector=unique(new_vector(:,1));
                nConditions=length(condition_vector);
                
                cond_matrix=zeros(nConditions,14);
                for iCond=1:nConditions
                    condition_nr=condition_vector(iCond);
                    info=MWorks_info(condition_nr);
                    frames=find(new_vector(:,1)==condition_nr)';
                    stim_data=[iCond condition_nr info.stim_present length(frames) frames([1 end]) info.condition_nr info.position_nr info.stim_id info.pos_x info.pos_y info.size_x info.size_y info.rotation];
                    cond_matrix(iCond,1:size(stim_data,2))=stim_data;
                end
                
                %%
                condition_vector=new_vector(:,1);
                nConditions=length(condition_vector);
                
                stim_matrix=zeros(nConditions,11);
                for iCond=1:nConditions
                    condition_nr=condition_vector(iCond);
                    info=MWorks_info(condition_nr);
                    frames=find(new_vector(:,1)==condition_nr)';
                    trial_nr=condition_nr-min(new_vector(:,1))+1;
                    if info.stim_present==1
                        stim_data=[iCond trial_nr info.stim_present info.condition_nr info.stim_id info.position_nr info.pos_x info.pos_y info.size_x info.size_y info.rotation];
                    else
                        stim_data=[iCond trial_nr info.stim_present ones(1,8)*-1];
                    end
                    stim_matrix(iCond,1:size(stim_data,2))=stim_data;
                end
                
                %% put new matrix in Experiment_info
                session_data.Experiment_info.cond_matrix=cond_matrix;
                if isfield(session_data.Experiment_info,'stim_matrix') && ~isfield(session_data.Experiment_info,'stim_matrix_bkp')
                    % store original in backup
                    session_data.Experiment_info.stim_matrix_bkp=session_data.Experiment_info.stim_matrix;
                end
                session_data.Experiment_info.stim_matrix=stim_matrix;
                
                
                if 0
                    %% sanity check for healthy file
                    A=session_data.Experiment_info.stim_matrix_bkp;
                    B=session_data.Experiment_info.stim_matrix;
                    all(eq(A(:),B(:)))
                end
                
                
                if are_you_sure==true
                    %% when we are sure this code is solid, save it to the data file
                    session_data.last_action='fix_missing_bitCodes';
                    session_data.updated=1;
                    session_data.save_data()
                else
                    %cond_matrix(:,[1 3:9])
                    session_data.Experiment_info
                end
            end
            
        end
    catch
        e=lasterror;
        disp(e.message)
    end
end



%% new thing
% we have the truth = MWorks, with stim properties and timestamp
% we can resample the bitcode to match the scim bitcodes,
% once match is found, we can couple frame to stim properties
% we will narrow it down to stim properties per condition and frame times
% start and end

%% to get there, we need 2 things:
% MW: a matrix that has bitcode, stim properties, event_id/condition_id
% SI: a matrix that has bitcode and frame_nr
% those two can be linked, so we can assign a condition_id to a set of
% frames

