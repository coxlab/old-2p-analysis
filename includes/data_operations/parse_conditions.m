function varargout=parse_conditions(varargin)

X=varargin{1};

%%% Prepend -1 to detect first condition
D=diff([-1;X]);
T_start=find(D>0);
T_end=find(D<0)-1;

apply_corrections=0;
if apply_corrections==1
    %%% Find and remove spurious on transitions
    x=find(D>0); % find negative transitions
    x(x<=1)=[];
    check=X(x-1); % and check previous
    T_start(check>0)=[]; % remove if these are positive
    
    %%% Find and remove spurious off transitions
    check=X(D<0); % find negative transitions
    T_end(check>0)=[]; % remove if these are positive
end

%%% Measure both vectors
N1=length(T_start);
N2=length(T_end);

%%% If equal, we can process
if N1==N2
    M=[(1:N1)' T_start T_end T_end-T_start+1 X(T_start) X(T_end) X(T_start)-X(T_end)];
    sel=diff(M(:,5:6),[],2)==0;
    M(sel==0,:)=[];
    if all(diff(M(:,5:6),[],2)==0)
        varargout{1}=M;
    else % show problem
        disp([M(:,5:6) diff(M(:,5:6),[],2)])
    end
else
    N=min([N1 N2]);
    M=[(1:N)' T_start(1:N) T_end(1:N) T_end(1:N)-T_start(1:N)+1 X(T_start(1:N)) X(T_end(1:N)) X(T_start(1:N))-X(T_end(1:N))];
    sel=diff(M(:,5:6),[],2)==0;
    M(sel==0,:)=[];
    if all(diff(M(:,5:6),[],2)==0)
        varargout{1}=M;
    else % show problem
        disp([M(:,5:6) diff(M(:,5:6),[],2)])
    end
    % else % if not: error
    %     N=min([N1 N2]);
    %     [T_start(1:N) T_end(1:N) T_end(1:N)-T_start(1:N)+1]
    %
    %     varargout{1}=-1;
end