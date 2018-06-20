function [ AT ] = AggregateNetwork(A)
% AggregateNetwork Aggregates multilayer network by summing over all layers

% Version: 2.0.2
% Date: Wed 20 Jun 2018 16:01:02 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com

T=length(A);
AT=A{1};
if T>1
    for i=2:T
        AT=AT+A{i};
    end
end

end

