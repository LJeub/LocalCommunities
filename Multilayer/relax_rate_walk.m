function [P,id,layer]=relax_rate_walk(A)
% relax_rate_walk compute relax-rate walk matrices
%
% Input: 
%
%   A: cell array of adjacency matrices for each layer
%
% Output: 
%
%   P: function that returns relax-rate-walk matrix given relax rate as
%      input
%
%   id: vector of node ids
%
%   layer: vector of layer ids   

% Version: 2.0
% Date: Mon 25 Jul 2016 17:06:57 BST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk
l_width=size(A{1},1);
n_l=length(A);
id=zeros(l_width*n_l,1);
layer=zeros(l_width*n_l,1);
pos=0;
for i=1:n_l
    [ii,jj,vv]=find(A{i});
    ASi(pos+(1:length(ii)))=ii+l_width*(i-1);
    ASj(pos+(1:length(jj)))=jj+l_width*(i-1);
    ARj(pos+(1:length(jj)))=jj;
    ASv(pos+(1:length(vv)))=vv;
    pos=pos+length(ii);
    id((1:l_width)+(i-1)*l_width)=1:l_width;
    layer((1:l_width)+(i-1)*l_width)=i;
end

AS=sparse(ASi,ASj,ASv,n_l*l_width,n_l*l_width);
AR=repmat(sparse(ASi,ARj,ASv,n_l*l_width,l_width),1,n_l);

ksout=sum(AS,1);
krout=sum(AR,1);
ksin=sum(AS,2);

PS=div_0(AS,repmat(ksout,l_width*n_l,1));
%for i=find(ksout(:)'==0)
%    PS(i,i)=1;
%end

PR=div_0(AR,repmat(krout,l_width*n_l,1));
%for i=find(krout(:)'==0)
%    PR(i,i)=1;
%end

    function P=walk_matrix(r)
        P=(1-r)*PS+r*PR;
    end
    
P=@walk_matrix;
end
