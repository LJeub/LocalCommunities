function [support, conductance, flag, connected]=ACLcut(W,d,seed,alpha,epsilon,max_vol)
% ACLcut local truncated PageRank method
%
% Implements the ACLcut truncated pagerank diffusion method to find 
% cuts with small conductance values around a seed vertex.
%   
% Inputs:
%   
%   W: adjacency matrix
%   d: vector of node strength
%   seed: vector of seed node indeces
%   alpha: value of teleportation parameter
%   epsilon: truncation parameter for APPR vector
%   max_vol: maximum volume of communities returned
%
% Outputs:
%
%   support: nodes in support of APPR vector in sweep cut ordering
%   conductance: vector of conductance values for each sweep set
%   flag: convergence flag (false if APPR vector calculation has not
%       converged after 10^6 iterations)
%   connected: identifies connected sweep sets
%
% Reference:
%
%   Andersen, R., Chung, F. R. K., & Lang, K. J. (2006). 
%   Local Graph Partitioning using PageRank Vectors . 
%   Proceedings of the 47th Annual Symposium on Foundations of Computer Science, 
%   IEEE (pp. 475?486). http://doi.org/10.1109/FOCS.2006.44
%
% see also NCP sweep_cut APPR

% Version: 2.0.2
% Date: Wed 20 Jun 2018 16:01:01 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com


% Up to half the volume of the network by default
if nargin<6
    max_vol=0.5*sum(d);
end

%unpack passed cell array 
if iscell(seed)
    if length(seed)==1
        seed=seed{1};
    else
        error('ACLcut only accepts a single seed set as input')
    end
end


%compute approximate pagerank vector for seed
r=zeros(size(d));
r(seed)=1/length(seed);
[p,flag]=APPR(alpha,epsilon,r,W,d);

%Compute normalized sweep vector
p(p>0)=p(p>0)./d(p>0);

%compute conductance for sweep sets
[conductance,support,connected]=sweep_cut(p,W,d,max_vol);

end

