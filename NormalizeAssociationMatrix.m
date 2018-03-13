function assoc_mat=NormalizeAssociationMatrix(assoc_mat)
% Normalize an association matrix returned by NCP
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

% Version: 2.0.1
% Date: Tue 13 Mar 2018 15:46:52 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

n=length(assoc_mat);
norm=repmat(diag(assoc_mat),1,n)+repmat(diag(assoc_mat)',n,1)-assoc_mat;

assoc_mat=assoc_mat./norm;
end
