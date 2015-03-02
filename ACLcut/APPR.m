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
%
% Implemented as a mex function

% Version: 1.01
% Date: Tue 25 Mar 2014 17:22:32 GMT
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

end