function varargout=smart_shuffle(varargin)

X=varargin{1};
nFrames=length(X);

%%% break up into stim and blank blocks
T_start=find([0;diff(X)>0]);
delete_last=0;
if X(end)>0
    T_end=[find(diff(X(2:end))<0)+1;nFrames];
    delete_last=1;
else
    T_end=find(diff(X(2:end))<0)+1;
end
condition_matrix=[T_start T_end X(T_start) X(T_end) T_end-T_start+1];

if 0
    %%% remove last trial when no blank upon exp end
    if delete_last==1
        condition_matrix(end,:)=[];
    end
end

if 0
    %%% remove condition that have more frames then requested
    sel=condition_matrix(:,5)==mode(condition_matrix(:,5));
    condition_matrix(sel==0,:)=[];
end

%%% sanity check whether begin and end match
check=all(eq(X(T_start),X(T_end)));
if check==0
    disp(condition_matrix)
    die
end

nTrials=size(condition_matrix,1);
n_rows=nFrames;
%n_rows=max(condition_matrix(:,2));
shuffled=condition_matrix(randperm(nTrials),3);

X_shuffled=zeros(n_rows,1)-1;
%X=zeros(n_rows,1)-1;
for iTrial=1:nTrials
    indices=condition_matrix(iTrial,1):condition_matrix(iTrial,2);
    %X(indices)=condition_matrix(iTrial,3);
    X_shuffled(indices)=shuffled(iTrial);
end

%[X X_shuffled]
%mean([X X_shuffled])
%n_rows

varargout{1}=X_shuffled;