function [ AT ] = AggregateNetwork(A)
% Aggregates multilayer network by summing over all layers

T=length(A);
AT=A{1};
if T>1
    for i=2:T
        AT=AT+A{i};
    end
end

end

