function [ AT ] = AggregateNetwork(A)
% AggregateNetwork Aggregates multilayer network by summing over all layers

% Version: 2.0
% Date: Mon 25 Jul 2016 17:06:57 BST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

T=length(A);
AT=A{1};
if T>1
    for i=2:T
        AT=AT+A{i};
    end
end

end

