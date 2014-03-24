function [conductance,support,connected]=sweep_cut(p,W,d,max_vol)
% [conductance,support,connected]=sweep_cut(p,W,d,max_vol)
% 
% computes a sweep cut for ranking vector p and adjacency matrix W with
% node strengths d. Optionally specify max_vol to restrict the maximum
% volume of sweep sets.
%
% The support is taken over positive values of p unless max_vol is
% specified, in which case support is taken over all nodes
%
% Outputs:
%
%   conductance: vector of conductance values for each sweep set
%   support: node indeces for sweep sets, such that the set of nodes given
%       by support(1:k) has conductance given by conductance(k)
%   connected: logical vector with elements indicating whether the subgraph
%       induced by the corresponding sweep set is connected

% Version:
% Date:
% Author:
% Email:

n=length(W);
mm=sum(d);

[ps,support]=sort(p,'Descend');

if nargin<4
%support
supp=find(ps<=0,1,'first')-1;
max_vol=inf;
else
supp=n-1;
end


if isempty(supp)
    supp=n-1;
end

if supp==0
    conductance=nan;
    support=[];
    connected=[];
else
    
    %logical indeces of nodes not in community
    indexcomp=true(n,1); 
    
    E=zeros(1,supp);
    E_p=0;
    i=1;
    Volsc=mm-d(support(1));
    Vols=d(support(1));
    conductance=ones(1,supp);
    connected=false(1,supp);
    
    components=zeros(size(W,1),1);
    
    num_comp=0;
    while i<=supp && Vols<=max_vol
   
        indexcomp(support(i))=false;
        
        %compute change in conductance by adding node index(i)
        E(i)=E_p-sum(sum(W(support(1:i-1),support(i))))+sum(sum(W(indexcomp,support(i))));
        E_p=E(i);
       
         %store conductance values
        conductance(i)=E(i)/min(Volsc,Vols);
		% cond(i)=mm*E(i)/(Volsc*Vols);
        
        % check if sweepset is connected
        if nargout>2
            n_comp=components(W(:,support(i))~=0);
            n_comp=n_comp(n_comp>0);
            
            if isempty(n_comp)
                components(support(i))=num_comp+1;
                num_comp=num_comp+1;
            else
                merge=unique(n_comp);
                components(support(i))=merge(1);
                if length(merge)>1
                    no_merge=true(num_comp,1);
                    no_merge(merge)=false;
                    ind=support(1:i);
                    C=sparse(1:i,components(ind),true);
                    Cn=[C(:,no_merge),any(C(:,merge),2)];
                    num_comp=size(Cn,2);
                    for it=1:num_comp
                        components(ind(Cn(:,it)))=it;
                    end
                end
            end
            connected(i)=num_comp==1;
        end
        
		%compute change in volumes
         Volsc=Volsc-d(support(i+1));
         Vols=Vols+d(support(i+1));
         i=i+1;
    end
    support=support(1:i-1);
    conductance=conductance(1:i-1);
    connected=connected(1:i-1);
end


end


