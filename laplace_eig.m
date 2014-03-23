function [lambda_2,V]=laplace_eig(A,tol,maxiter)
% compute second-smallest eigenvalue of normalised Laplacian matrix
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
    
    R=ichol(L);
    RT=R';
    precfun = @(x)RT\(R\x);
    
    in2=rand(length(A),1);
    [V,lambda_2,flag]=lobpcg(in2,L,[],precfun,tol,maxiter,0,d);
    if flag
        warning('not converged')
    end
end