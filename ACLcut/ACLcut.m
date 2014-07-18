function [support, conductance, flag, connected]=ACLcut(W,d,seed,alpha,epsilon)
% [support,conductance,volume]=ACLcut(W,d,seed,alpha,epsilon)
%
% ACLcut: implements the ACLcut truncated pagerank diffusion method to find small
%   conductance cuts around a seed vertex
%   
% Inputs:
%   
%   W: adjacency matrix
%   d: vector of node strength
%   seed: seed node
%   alpha: value of teleportation probability
%   epsilon: truncation parameter for APPR vector
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

% Version: 1.02
% Date: Fri 18 Jul 2014 13:16:55 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

%compute approximate pagerank vector for seed
[p,flag]=APPR(alpha,epsilon,seed,W,d);

%Compute normalized sweep vector
p=p(:)./d(:);

%compute conductance for sweep sets
[conductance,support,connected]=sweep_cut(p,W,d);
if isempty(support)
    support=seed;
    conductance=1;
    connected=true;
end

end

