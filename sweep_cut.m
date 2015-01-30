function [conductance,support,connected]=sweep_cut(p,W,d,max_vol)
% [conductance,support,connected]=sweep_cut(p,W,d,max_vol)
% 
% sweep_cut: computes a sweep cut for ranking vector p and adjacency matrix W with
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
%
% see also ACLcut MOVcut EGOcut global_cut

% Version: 1.01
% Date: Tue 25 Mar 2014 17:22:33 GMT
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

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

supp=min(supp,ceil(n/2));


if supp==0
    conductance=nan;
    support=[];
    connected=[];
else
    
    %logical indeces of nodes not in community
    %indexcomp=true(n,1); 
    
    E=0;
    i=1;
    Volsc=mm-d(support(1));
    Vols=d(support(1));
    conductance=ones(1,supp);
    connected=false(1,supp);
    
    if nargout>2
        components=zeros(supp,1);
        WS=W(support(1:supp),support(1:supp));
    end
    
    num_comp=0;
    while i<=supp && Vols<=max_vol
   
        %indexcomp(support(i))=false;
        
        %compute change in conductance by adding node index(i)
        E=E-sum(sum(W(support(i),support(1:i-1))))+sum(sum(W(support(i+1:end),support(i))));

       
         %store conductance values
        conductance(i)=E/min(Volsc,Vols);
		% cond(i)=mm*E(i)/(Volsc*Vols);
        
        % check if sweepset is connected
        if nargout>2
            merge=unique_fast([components(WS(1:i-1,i)~=0);components(WS(i,1:i-1)~=0)]);
           
            if isempty(merge)
                components(i)=num_comp+1;
                num_comp=num_comp+1;
            else
                %merge=unique(merge);
                components(i)=merge(1);
                if length(merge)>1
                    no_merge=true(num_comp,1);
                    no_merge(merge)=false;
                    ind=1:i;
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


