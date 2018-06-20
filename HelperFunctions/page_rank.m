function p=page_rank(P,alpha,s)
% Compute page-rank vector
%
% Input:
%           P: transition matrix
%
%           alpha: teleportation parameter (if alpha=0 computes stationary
%                  distribution of P
%
%           s: seed vector (defaults to uniform distribution)
%
% Output:
%           p: page-rank vector

if nargin<3
    s=ones(length(P),1)/length(P);
else
    s=s/sum(s);
end

if alpha==0
    [p,~]=eigs(P,1,'la');
else
    p=bicgstab((speye(length(P))+(alpha-1)*P),alpha*s,[],100);
end

% correct for potential numerical error (make sure p satisfies assumptions)
if sum(p)<0
    p=-p;
end
p=max(0,p);
p=p/sum(p);

end
