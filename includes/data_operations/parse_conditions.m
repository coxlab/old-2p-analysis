function varargout=parse_conditions(varargin)

X=varargin{1};

D=diff([-1;X]);

T_start=find(D>0);

T_end=find(D<0)-1;
check=X(D<0);
T_end(check>0)=[];

N1=length(T_start);
N2=length(T_end);

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
    %[T_start(1:N) T_end(1:N) T_end(1:N)-T_start(1:N)+1]
    
    varargout{1}=-1;
end