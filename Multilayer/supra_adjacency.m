function [AS,id,layer]=supra_adjacency(A,omega)
% supra_adjacency Convert multilayer network to supra-adjacency matrix
%
% Input:
%
%   A: cell array of adjacency matrices for each layer
%
%   omega: weight of interlayer edges
%
% Output:
%
%   AS: supra-adjacency matrix with interlayer edges of weight omega with
%       all-to-all (multiplex) coupling
%
%   id: vector of node ids for state nodes that are present (i.e. have
%       connections)
%
%   layer: vector of layer ids for state nodes that are present

% Version: 2.0.2
% Date: Wed 20 Jun 2018 16:01:02 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com
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
