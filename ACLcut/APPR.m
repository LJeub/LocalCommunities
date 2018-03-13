function [p,not_converged,r]=APPR(alpha,epsilon,s,A,d)
% Compute a personalized PageRank vector
% 
% Implements the ApproximatePR algorithm from:
%   Andersen, R., Chung, F. R. K., & Lang, K. J. (2006).
%   Local Graph Partitioning using PageRank Vectors (pp. 475-486).
%   FOCS'06
%
% Input:
% 
%   alpha: teleportation parameter between 0 and 1
%   epsilon: truncation parameter
%   s: seed node
%   A: adjacency matrix
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

% Version: 2.0.1
% Date: Tue 13 Mar 2018 15:46:51 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

% This function should never run if everything is installed correctly
error(['APPR mex-function not found, make sure to compile the C++ code by '...
    'running the ''compile_APPR_mex'' script'])
end

