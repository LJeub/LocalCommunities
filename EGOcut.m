function [support,conductance,flag,sweep_set]=EGOcut(W,d,seed,~,~)
% [support,conductance,flag,sweep_set]=EGOcut(W,d,seed,~,~)
%
% EGOcut: computes sweep cuts based on ranking nodes by geodesic 
%   distance from a seed node.
%
% Inputs: 
%
%   W: adjacency matrix
%   d: vector of node strength
%   seed: seed node
%
% Outputs:
%   
%   support: node indeces in sweep cut ordering
%   conductance: vector of conductance values for each sweep set
%   flag: always false (no convergence issues to deal with)
%   sweep_set: identifies sets that correspond to sweep sets, i.e.,
%       where all nodes with egorank equal to the minimum in the 
%       community have been included.

% Version: 1.01
% Date: Tue 25 Mar 2014 17:22:32 GMT
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

% use inverse weight as distance
D=W;
D(W>0)=1./(W(W>0));
egorank=shortest_paths(D,seed);
egorank=1./(1+egorank);

[conductance,support]=sweep_cut(egorank,W,d);

%sweep_set indicates actual communities, where all nodes
%with egorank equal to the minimum in the community have been included.
%(All communiites returned by this procedure are connected)
sweep_set=logical(diff(sort(egorank,'Descend')));
sweep_set=sweep_set(:)';


% no convergence issues to indicate
flag = false;

end


