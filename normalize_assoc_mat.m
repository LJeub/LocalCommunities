function assoc_mat=normalize_assoc_mat(assoc_mat)
% [assoc_mat]=normalize_assoc_mat(assoc_mat)
%
% normalize_assoc_mat: normalize an association matrix returned by NCP
% 
% Normalize an association matrix assoc_mat returned by NCP by dividing 
% each entry by the total number of times either node has appeared in a
% community
%
% Input: 
%   assoc_mat: association matrix (last output from NCP)
%
% Output:
%   assoc_mat: normalized version of the association matrix
%
% see also NCP

% Version: 2.0-beta
% Date: Fri 17 Jun 2016 17:33:45 BST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

n=length(assoc_mat);
norm=repmat(diag(assoc_mat),1,n)+repmat(diag(assoc_mat)',n,1)-assoc_mat;

assoc_mat=assoc_mat./norm;
end
