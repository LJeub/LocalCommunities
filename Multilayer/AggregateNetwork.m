function [ AT ] = AggregateNetwork(A)
% AggregateNetwork Aggregates multilayer network by summing over all layers

% Version: 2.0.1
% Date: Tue 13 Mar 2018 15:46:52 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

T=length(A);
AT=A{1};
if T>1
    for i=2:T
        AT=AT+A{i};
    end
end

end

