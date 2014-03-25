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

% Version: 1.0
% Date: Tue 25 Mar 2014 16:10:49 GMT
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

[~,v]=laplace_eig(A);

[cond,ind]=sweep_cut(v,A,sum(A,2),inf);

[cond,k]=min(cond);

S=ones(length(A),1)*2;

S(ind(1:k))=1;
end
