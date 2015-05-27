function varargout=analyse_RF(varargin)
Y=varargin{1};
X=varargin{2};

%%% use the fact that blanks occur much more frequently to exclude these
%%% from analysis
condition_vector=unique(X);
counts=hist(X,length(condition_vector));
condition_vector=condition_vector(between(zscore(counts),[-1 1]));
nConditions=length(condition_vector);
nFrames=size(X,1);

%%% Safe frames before end of experiment
stim_length_frames=24;

%%% Which frames relative to stim onset do we use
frame_selector=[2 7];

%%% Loop over conditions
for iCond=1:nConditions
    condition_nr=condition_vector(iCond);
    
    %%% find different presentation of condition
    frames_selected=find(X==condition_nr);
    start_frames=[1;diff(frames_selected)>1];
    repeat_starts=frames_selected(start_frames==1);
    repeat_starts(repeat_starts+stim_length_frames>nFrames)=[];
    nRepeats=length(repeat_starts);
    
    %%% Collect mean response for each repeat
    repeat_vector=zeros(nRepeats,1);
    for iRepeat=1:nRepeats
        repeat_start=repeat_starts(iRepeat);
        repeat_vector(iRepeat)=mean(Y(repeat_start+frame_selector(1):repeat_start+frame_selector(2)));
        %stim_matrix(repeat_start-2:repeat_start+6,:); % sanity check
    end
    
    result.nRepeats=nRepeats;
    result.mean_resp=mean(repeat_vector);
    result.std_resp=std(repeat_vector);
    
    
end

varargout{1}=result;