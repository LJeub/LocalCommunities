function [AS,id,layer]=supra_adjacency(A,omega)
N=size(A{1},1);
T=length(A);

pos=0;
for i=1:length(A)
    cid=find(sum(A{i},1)|sum(A{i},2)');
    id(pos+(1:length(cid)))=cid;
    layer(pos+(1:length(cid)))=i;
    AS((i-1)*N+(1:length(A{i})),(i-1)*N+(1:length(A{i})))=A{i};
    pos=pos+length(cid);
end

G=sparse(id+(layer-1)*N,id,true);
for i=1:size(G,2)
    AS(G(:,i),G(:,i))=omega;
end

AS=AS-diag(diag(AS));
end
