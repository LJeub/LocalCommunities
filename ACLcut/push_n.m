function [ p,r] = push_n( p,r,u,alpha,W,d,n )
%PUSH_N apply push n times at once

% Version: 0.1-beta
% Date: Mon 24 Mar 2014 21:39:53 GMT
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

mult=(1-((1-alpha)/2).^n)/(1-(1-alpha)/2);
p(u) = p(u) + alpha*r(u)*mult;
r=r + (1-alpha)*r(u)/(2*d(u))*mult*W(:,u);
r(u) = (1-alpha)^n*r(u)/(2^n);

end

