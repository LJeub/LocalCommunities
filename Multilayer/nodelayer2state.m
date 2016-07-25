function state=nodelayer2state(N,nodelayer)
% convert node-layer indeces to state indeces
%
% Input:
%           N: number of physical nodes of the network
%
%           nodelayer: nx2 matrix of node-layer indeces or cell array of
%                      nx2 matrices of node-layer indeces.
%
% Output:
%           state: nx1 vector of state indeces
%
% Note: if input is a vector, the function returns the input unchanged.

% Version: 2.0
% Date: Mon 25 Jul 2016 17:06:57 BST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

if iscell(nodelayer)
    state=cell(size(nodelayer));
    for nl=1:numel(nodelayer)
        state{nl}=nodelayer2state(N,nodelayer{nl});
    end
else
    if size(nodelayer,2)==2
        state=(nodelayer(:,2)-1)*N+nodelayer(:,1);
    elseif isvector(nodelayer)
        state=nodelayer;
    else
        error('Format of nodelayer input is invalid');
    end
end

end
