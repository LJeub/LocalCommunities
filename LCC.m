function [ A_c,ind] = LCC(A )
%LCC [A_c,ind]=LCC(A) find the largest connected component of a network
%
%   Input:
%       A: adjacency matrix
%
%   Outputs:
%       A_c: adjacency matrix of largest connected component
%       ind: node indeces of nodes in the largest connected component, i.e.
%           A(ind,ind)=A_c

% Lucas Jeub
% jeub@maths.ox.ac.uk

[C,sizes]=components(A);

[~,k]=max(sizes);
ind=find(C==k);

A_c=A(ind,ind);

end

