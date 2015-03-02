function [support, conductance, flag, connected]=ACLcut(W,d,seed,alpha,epsilon,max_vol)
% [support,conductance,volume]=ACLcut(W,d,seed,alpha,epsilon)
%
% ACLcut: implements the ACLcut truncated pagerank diffusion method to find small
%   conductance cuts around a seed vertex
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
% see also APPR sweep_cut NCP MOVcut

% Version: 1.01
% Date: Tue 25 Mar 2014 17:22:32 GMT
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk


% no volume restriction by default
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
%[p,flag]=APPR_nomex(alpha,epsilon,seed,W,d);
r=zeros(size(d));
r(seed)=1/length(seed);
[p,flag]=APPR(alpha,epsilon,r,W,d);

%Compute normalized sweep vector
p=div_0(p(:),d(:));

%compute conductance for sweep sets
[conductance,support,connected]=sweep_cut(p,W,d,max_vol);


end

