function [mapping, labels]=getMapping(string_array)

labels=unique(string_array);
N=length(labels);
mapping=zeros(N,1);
for index=1:N
    sel=ismember(string_array,labels{index});
    mapping(sel)=index;
end
mapping=mapping(:);
