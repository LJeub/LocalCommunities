function nodelayer=state2nodelayer(N,state)
% convert state indeces to node-layer indeces
%
% Input:
%           N: number of nodes in the network
%
%           state: vector of state indeces or cell array of
%                  vectors of state indeces.
%
% Output:
%           nodelayer: matrix of node-layer indeces

% Version: 2.0-beta
% Date: Fri 17 Jun 2016 17:33:45 BST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

if iscell(state)
    nodelayer=cell(size(state));
    for st=1:numel(state)
        nodelayer{st}=state2nodelayer(N,state{st});
    end
else
    nodelayer=zeros(numel(state),2);
    nodelayer(:,1)=mod(state(:)-1,N)+1;
    nodelayer(:,2)=floor((state(:)-1)./N)+1;
end

end
