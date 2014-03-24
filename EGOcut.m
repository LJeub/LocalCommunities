function [support, conductance, flag, connected]=EGOcut(W,d,seed,~,~)
%

% Version:
% Date:
% Author:
% Email:

% use inverse weight as distance
D=W;
D(W>0)=1./(W(W>0));
egorank=shortest_paths(D,seed);
egorank=1./(1+egorank);

[conductance,support]=sweep_cut(egorank,W,d);

%this is a hack, connected indicates actual communities, where all nodes
%with egorank equal to the minimum in the community have been included.
%(All communiites returned by this procedure are connected)
connected=logical(diff(sort(egorank,'Descend')));
connected=connected(:)';


% no convergence issues to indicate
flag = false;

end


