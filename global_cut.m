function [S,cond]=global_cut(A)

% Version: 0.1-beta
% Date: Mon 24 Mar 2014 21:39:53 GMT
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

[~,v]=laplace_eig(A);

[cond,ind]=sweep_cut(v,A,sum(A,2),inf);

[cond,k]=min(cond);

S=ones(length(A),1)*2;

S(ind(1:k))=1;
end
