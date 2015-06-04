function response=between(values,limits,varargin)
%function response=between(values,limits,varargin)
% check whether VALUES fall between interval defined by LIMITS
% works for 2 column matrix with start and stop values of the interval.

if nargin==0
    help between
elseif nargin==1
    response=values>=0&values<1;
elseif nargin==2
    % serialize values
    valuesSize=size(values);
    values=reshape(values,1,numel(values));
    response=zeros([valuesSize size(limits,1),]);
    for limit_index=1:size(limits,1)
        responses=values>=limits(limit_index,1)&values<limits(limit_index,2);
        % reshape to plane if nec
        response_temp=reshape(responses,valuesSize(1),valuesSize(2));
        response(:,:,limit_index)=response_temp;
    end
    response=logical(response);
elseif nargin==3
    borders=varargin{1};
    L=limits(1);
    U=limits(2);
    
    if eq(borders,[0 0])
        response=values>L&values<U;
    elseif eq(borders,[0 1])
        response=values>L&values<=U;
    elseif eq(borders,[1 0])
        response=values>=L&values<U;
    elseif eq(borders,[1 1])
        response=values>=L&values<=U;
    end
end

