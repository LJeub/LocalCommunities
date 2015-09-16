function [S,cond]=global_cut(A)
%[S,cond]=global_cut(A)
%
% global_cut: approximate the globally minimal conductance
%   cut
%
% Input:
%   A: adjacency matrix
%
% Outputs:
%   S: indicates bipartition found (entries are 1 or 2,
%       indicating which part of the partition the corresponding
%       node belongs to)
%   cond: the conductance value of the bipartition

% Version: 1.02
% Date: Fri 18 Jul 2014 13:16:55 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

[~,v]=laplace_eig(A);

[cond,ind]=sweep_cut(v,A,sum(A,2),inf);

[cond,k]=min(cond);

S=ones(length(A),1)*2;

S(ind(1:k))=1;
end
