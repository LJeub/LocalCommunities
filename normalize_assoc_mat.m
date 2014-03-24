function assoc_mat=normalize_assoc_mat(assoc_mat)
% NORMALIZE_ASSOC_MAT normalize an association matrix returned by NCP
% 
% [assoc_mat] = normalize_assoc_mat(assoc_mat)
% 
% Normalize an association matrix assoc_mat returned by NCP by dividing 
% each entry by the total number of times either node has appeared in a
% community
%
% see also NCP

% Version: 0.1-beta
% Date: Mon 24 Mar 2014 21:39:53 GMT
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

n=length(assoc_mat);
norm=repmat(diag(assoc_mat),1,n)+repmat(diag(assoc_mat)',n,1)-assoc_mat;

assoc_mat=assoc_mat./norm;
end
