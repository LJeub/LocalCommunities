function [support, conductance, flag, connected]=MOVcut(W,d,seed,gamma,c)
% [support, conductance, flag, connected]=MOVcut(W,d,seed,gamma,c)
%
% MOVcut: implements the MOVcut locally biased spectral optimisation to find low
% conductance cuts around a seed vertex, see:
%   Mahoney, M. W., Orecchia, L., & Vishnoi, N. K. (2012). 
%   A local spectral method for graphs: With applications to improving graph
%       partitions and exploring data graphs locally. JMLR, 13, 2339?2365.
%
% Input:
%   
%   W: adjacency matrix
%   d: vector of node strength
%   seed: nodes of interest
%   gamma: value of teleportation parameter
%   c: volume factor (if c is not specified, support is only over
%       positive values of the GPPR vector)
%
% Output:
%
%   support: nodes  in sweep cut ordering
%   conductance: vector of conductance values for each sweep set
%   flag: convergence flag
%   connected: identifies connected sweep sets
%
% see also GPPR sweep_cut NCP ACLcut

% Version: 0.1-beta
% Date: Mon 24 Mar 2014 21:39:53 GMT
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

D=diag(d);

s=sparse(seed,1,1,size(W,2),1); %creates seed vector using seed
s=s-(s'*d(:)/(d(:)'*d(:)))*d(:); %othogonalizes s relative to degree sequence
s=s/sqrt((s'*D*s)); %normalizes s appropriately

[p,flag]=GPPR(gamma,s,W,d); %Computes the GPPR vector

if nargin>4
    
    p=p/sqrt(p'*D*p);
    kappa=(s'*D*p)^2;
    max_vol=c/kappa;
    
    [conductance,support,connected]=sweep_cut(p,W,d,max_vol);
    
else
    [conductance,support,connected]=sweep_cut(p,W,d);
end
