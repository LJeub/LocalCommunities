function [support,conductance,flag,sweep_set]=EGOcut(W,d,seed,~,~,max_vol)
% Compute sweep cuts based on ranking nodes by geodesic distance from a seed node.
%
% Inputs: 
%
%   W: adjacency matrix
%   d: vector of node strength
%   seed: seed node or seed set of nodes
%
% Outputs:
%   
%   support: node indeces in sweep cut ordering
%   conductance: vector of conductance values for each sweep set
%   flag: always false (no convergence issues to deal with)
%   sweep_set: identifies sets that correspond to sweep sets, i.e.,
%       where all nodes with egorank equal to the minimum in the 
%       community have been included.
%
% see also NCP sweep_cut

% Version: 2.0.2
% Date: Wed 20 Jun 2018 16:01:02 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com

% up to half the volume of the network by default
if nargin<6
    max_vol=0.5*sum(d);
end

if iscell(seed)
    if length(seed)==1
        seed=seed{1};
    else
        error('EGOcut only accepts a single seed set as input')
    end
end

% use inverse weight as distance
D=W;
D(W>0)=1./(W(W>0));
egorank=zeros(length(seed),length(W));
for i=1:length(seed)
    egorank(i,:)=shortest_paths(D,seed(i));
end
egorank=min(egorank,[],1);
egorank=1./(1+egorank);

[conductance,support]=sweep_cut(egorank,W,d,max_vol);

%sweep_set indicates actual communities, where all nodes
%with egorank equal to the minimum in the community have been included.
%(All communiites returned by this procedure are connected)
sweep_set=logical(diff(sort(egorank,'Descend')));
sweep_set=sweep_set(1:length(support));
sweep_set=sweep_set(:)';

% no convergence issues to indicate
flag = false;

end


