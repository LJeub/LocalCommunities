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

% Version: 0.1-beta
% Date: Mon 24 Mar 2014 21:39:53 GMT
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

[C,sizes]=components(A);

[~,k]=max(sizes);
ind=find(C==k);

A_c=A(ind,ind);

end

