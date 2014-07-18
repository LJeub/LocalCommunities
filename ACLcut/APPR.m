function [p,not_converged,r]=APPR(alpha,epsilon,v,W,d)
% [p,not_converged,r]=APPR(alpha,s,A)
%
% APPR: computes the personalized pagerank vector
% 
% Implements the ApproximatePR algrotihm from:
%   Andersen, R., Chung, F. R. K., & Lang, K. J. (2006).
%   Local Graph Partitioning using PageRank Vectors (pp. 475?486).
%   FOCS'06
%
% Input:
% 
%   alpha: teleportation parameter between 0 and 1
%   epsilon: truncation parameter
%   v: seed node
%   W: adjacency matrix
%   d: vector of node strengths
%
% Output:
%
%   p: PageRank vector as a row vector
%   not_converged: flag indicating that maxiter has been reached
%   r: residual vector
%
% see also ACLcut

% Version: 1.02
% Date: Fri 18 Jul 2014 13:16:55 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

maxiter=10^6;
not_converged=false;
n=length(d);

%creates seed vector with vertex seed v
r=zeros(n,1);
r(v)=1;

p=zeros(1,n);
resvec=r./d(:);

iter=0;
I=find(resvec>=epsilon);
while ~isempty(I)&&iter<maxiter
    k=length(I);
    iter =iter+k;
    for c=1:k
        n_iter=ceil(log(r(I(c))/(epsilon*d(I(c))))/log(2/(1-alpha))+eps);
        [p,r]=push_n(p,r,I(c),alpha,W,d,n_iter);
    end
    resvec=r./d(:);
    I=find(resvec>=epsilon);
    [~,s]=sort(resvec(I));
    I=I(s);
end

if iter>=maxiter
    warning('not converged after %f iterations, max residual still %f', iter,max(resvec))
    not_converged=true;
end
p=sparse(p);
end
