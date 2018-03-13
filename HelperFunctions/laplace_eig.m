function [lambda_2,V]=laplace_eig(A,tol,maxiter)
% Compute second-smallest eigenvalue and corresponding eigenvector of normalized Laplacian matrix
%
% Inputs: 
%   A: adjacency matrix
%   tol: error tollerance for the computation
%   maxiter: maximum number of iterations for the computation
%
% Outputs:
%   lambda_2: second-smallest eigenvalue of the normalized 
%       Laplacian matrix of A
%   V: corresponding eigenvector
%
% This function relies on the lobpcg method:
%   A. V. Knyazev, Toward the Optimal Preconditioned Eigensolver:
%   Locally Optimal Block Preconditioned Conjugate Gradient Method,
%   SIAM Journal on Scientific Computing 23 (2001), no. 2, pp. 517-541. 
%   http://dx.doi.org/10.1137/S1064827500366124
%
% see also lobpcg

% Version: 2.0.1
% Date: Tue 13 Mar 2018 15:46:51 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

if nargin<2
    tol=10^-6;
end

if nargin<3
    maxiter=1000;
end

if ~isequal(A,A')
    A=(A+A')/2;
    warning('symmetrised adjacency matrix')
end


D=diag(sum(A));
d=full(diag(D)).^(0.5);
L=(D-A);
D2=diag(diag(D).^-0.5);
L=D2*L*D2;
L=(L+L')/2; % eliminate roundoff error
if length(L)<50
    L=full(L);
    [V,lambda_2]=eig(L);
    lambda_2=diag(lambda_2);
    [lambda_2,s]=sort(lambda_2);
    lambda_2=lambda_2(2);
    V=V(:,s(2));
else
    try
        R=ichol(L);
        RT=R';
        precfun = @(x)RT\(R\x);
    catch er
        warning('Incomplete Cholesky factorisation failed')
        precfun = @(x) x;
    end
    
    in2=rand(length(A),1);
    [V,lambda_2,flag]=lobpcg(in2,L,[],precfun,tol,maxiter,0,d);
    if flag
        warning('not converged')
    end
end
