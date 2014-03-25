function [p,flag]=GPPR(gamma,s,W,d)
% [p,flag]=GPPR(gamma,s,W)
%
% GPPR: compute a Generalized Personalized Page Rank vector
%
% Inputs:
%
%   gamma: a number between (-infty, lambda_2(G)) where lambda_2(G) is the
%       second smallest eigenvalue of the normalized Laplacian
%   s: a seed vector
%   W: Adjacency matrix for G
%   d: vector of node strengths
%
% Outputs:
%
%   p: A generalized personal pagerank vector
%   flag: flag indicating convergence of bicgstab

% Version: 1.0
% Date: Tue 25 Mar 2014 16:10:49 GMT
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

D=sparse(diag(d)); %Construct diagonal matrix of degree sequence
Lcomb=D-W; %Constructs Combinatorial Laplacian

try
M=ichol(Lcomb-gamma*D); %incomplete Choleski factorization as preconditioner
catch 
    M=sparse(diag(diag(Lcomb-gamma*D)));
    disp('Incomplete Cholesky failed, using diagonal instead')
end
[p,flag,relres,iter]=bicgstab(Lcomb-gamma*D,D*s,10^-5,100,M,M');

%display shorter error messages
pcgmessages={'maxiter reached','ill-conditioned','stagnated','too small or large'};
if flag>0
    disp(['pcgerror = ',pcgmessages{flag}])
    disp(['gamma = ',num2str(gamma)])
    disp(['relres = ',num2str(relres)])
    disp(['iter = ',num2str(iter)])
end
