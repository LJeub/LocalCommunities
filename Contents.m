% LOCALCOMMUNITIES
%
% This code implements the local community detection methods used in
% 
% 1. Jeub, L. G. S., Balachandran, P., Porter, M. A., Mucha, P. J., & Mahoney, M. W. (2015). 
%    Think locally, act locally: Detection of small, medium-sized, and large communities in large networks. 
%    Physical Review E, 91(1), 012821. 
%    http://doi.org/10.1103/PhysRevE.91.012821
%
% 2. Jeub, L. G. S., Mahoney, M. W., Mucha, P. J., & Porter, M. A. (2015). 
%    A Local Perspective on Community Structure in Multilayer Networks. 
%    arXiv:1510.05185 [cs.SI].
%    http://arxiv.org/abs/1510.05185
%
% Main Interface
%
%   NCP             - Compute the NCP (Network Community Profile) for a network
%   MultilayerNCP   - Convenience wrapper around NCP for multiplex networks
%
% Further Diagnostics
%
%   InternalConductance         - Compute internal conductance of communities
%   NormalizeAssociationMatrix  - Normalize an association matrix returned by NCP
%
% See also ACLcut EGOcut MOVcut HelperFunctions Multilayer
