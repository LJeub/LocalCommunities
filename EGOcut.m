function [support, conductance, flag, connected]=EGOcut(W,d,seed,~,~)
%

% Version: 0.1-beta
% Date: Mon 24 Mar 2014 21:39:53 GMT
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

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


