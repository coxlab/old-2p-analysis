function varargout=analyse_RF(varargin)
Y=varargin{1};
X=varargin{2};

if length(X)~=length(Y)
    error('Both vectors should be equal length...')
end

%%% use the fact that blanks occur much more frequently to exclude these
%%% from analysis
condition_vector=unique(X);
counts=hist(X,length(condition_vector));
condition_vector=condition_vector(between(zscore(counts),[-1 1]));
nConditions=length(condition_vector);
if nConditions<32
    nConditions=32;
    %fprintf('Not all conditions were presented: %d/%d\n',[length(condition_vector) nConditions])
    condition_vector=1:32;    
end
nFrames=size(X,1);

%%% Safe frames before end of experiment
stim_length_frames=24;

%%% Which frames relative to stim onset do we use
frame_selector=[2 7];

%%% Loop over conditions
result=struct;
for iCond=1:nConditions
    condition_nr=condition_vector(iCond);
    
    %%% find different presentation of condition
    frames_selected=find(X==condition_nr);
    if sum(frames_selected)==0
        result(iCond).nRepeats=0;
        result(iCond).mean_resp=NaN;
        result(iCond).std_resp=NaN;
    else
        start_frames=[1;diff(frames_selected)>1];
        repeat_starts=frames_selected(start_frames==1);
        repeat_starts(repeat_starts+stim_length_frames>nFrames)=[];
        nRepeats=length(repeat_starts);
        
        if nRepeats>0
            %%% Collect mean response for each repeat
            repeat_vector=zeros(nRepeats,1);
            for iRepeat=1:nRepeats
                repeat_start=repeat_starts(iRepeat);
                mean_resp=mean(Y(repeat_start+frame_selector(1):repeat_start+frame_selector(2)));
                if isnan(mean_resp)
                    die
                end
                repeat_vector(iRepeat)=mean_resp;
                %stim_matrix(repeat_start-2:repeat_start+6,:); % sanity check
            end
            
            result(iCond).nRepeats=nRepeats;
            result(iCond).mean_resp=mean(repeat_vector);
            result(iCond).std_resp=std(repeat_vector);
        end
    end
end

varargout{1}=result;