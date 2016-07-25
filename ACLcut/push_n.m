function [ p,r] = push_n( p,r,u,alpha,W,d,n )
%PUSH_N apply push n times at once

% Version: 2.0-beta
% Date: Fri 17 Jun 2016 17:33:45 BST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

mult=(1-((1-alpha)/2).^n)/(1-(1-alpha)/2);
p(u) = p(u) + alpha*r(u)*mult;
r=r + (1-alpha)*r(u)/(2*d(u))*mult*W(:,u);
r(u) = (1-alpha)^n*r(u)/(2^n);

end

