function [conductance_con,communities_con,conductance_dis,communities_dis,assoc_mat]=NCP(W,cut_function,varargin)
% Compute the NCP (Network Community Profile) for a network
%
% This function (with the appropriate options listed below) can be used
% to compute approximate local and global NCPs for a network using any
% of the three methods in:
%
%   Jeub, L. G. S., Balachandran, P., Porter, M. A., Mucha, P. J.,
%   & Mahoney, M. W. (2014).
%   Think Locally, Act Locally: The Detection of Small, Medium-Sized, and
%   Large Communities in Large Networks. arXiv:1403.3795 [cs.SI]
%
% Inputs:
%           W: adjacency matrix
%
%           cut_function: choose method to find local communities
%               (either 'ACL', 'MOV', or 'EGO').
%
%           Additional options can be given using 'key',value pairs listed
%           below (with default values given in []).
%
% Outputs:
%           conductance_con: vector of minimum conductance values
%               considering only connected communities, e.g.
%               conductance_con(10) gives the minimum conductance found for
%               communities with 10 nodes that are connected.
%               For cut_function='EGO', this only considers communities
%               that are actual sweep sets (i.e. nodes with the same
%               distance from the seed node are either all included or not
%               included)
%
%           communities_con: returns communities that achieve the
%               minimum conductance values in conductance_con
%
%           conductance_dis: vector of minimum conductance values, also
%               allowing disconnected communities, and in the case 'EGO',
%               also communities that are not sweep sets.
%
%           communities_dis: communities corresponding to conductance_dis
%
%           assoc_mat: association matrix, assoc_mat(i,j) = number of times
%               nodes i and j have appeared together in a sampled
%               community. (The association matrix counts only the best
%               community for a seed node and choice of parameter values)
%
%
% Options:
%           nodes [all nodes]: number of nodes to sample for each pair of
%               parameter values (nodes are sampled uniformly at random).
%
%           local [[]]: compute a local NCP by specifying a vector of node
%               indeces (i.e. compute a NCP where only the specified nodes
%               are used as seed nodes). If the number of nodes to sample
%               given by the nodes option is less than the number of nodes
%               specified by local, nodes are sampled uniformly at random
%               from those specified in local.
%
%           alpha [network and method dependend]: teleportation parameter,
%               for 'ACL' this corresponds to the teleportation parameter
%               in the lazy random walk, and for 'MOV' this corresponds to
%               'gamma'. This option has no effect for 'EGO'.
%
%
%           truncation [network and method dependend]: for 'ACL' this is
%               the 'epsilon' parameter which controls the truncation in
%               the approximation, whereas for 'MOV' this is the 'c'
%               parameter, which limits the volume of communities. This
%               option has no effect for 'EGO'.
%
%           viscount [10]: minimum number of times each node needs to be in
%               the best community before the sampling is stopped.
%
%           aggressive: [true]: if set to `true`, a node that has
%               been in the best community at least `viscount` times is
%               never used as a seed node, if set to 'false', only stop 
%               iterating early when all nodes have been visited at least 
%               'viscount' times.
%
%           teleportation: [0.1]: for directed networks use unrecorded link
%               teleportations with teleportation parameter `teleportation`
%               to estimate stationary distribution (No effect when
%               `stationarydistribution` is given)
%
%           stationarydistribution: []: specify to use stationary
%               distribution to estimate adjacency matrix based on
%               unrecorded teleportations
%
%           transitionmatrix: [false]: if set to true, `A` is treated as a
%               random walk transition matrix instead of an adjacency
%               matrix
%
% For example, to compute an NCP for a network with adjacency matrix W
% using the ACLcut method with alpha=0.01, one would call the function as
%   conductance=NCP(W,'ACL','alpha',0.01);
%
% One can then plot the NCP using
%   loglog(conductance)
%
% see also ACLcut EGOcut MOVcut NormalizeAssociationMatrix

% Version: 2.0-beta
% Date: Fri 17 Jun 2016 17:33:45 BST
% Author: Lucas Jeub
% Email: jeub@maths.ox.ac.uk

% Parse Options
options=OptionStruct('nodes',length(W),'local',[],'alpha',[],...
    'truncation',[],'viscount',10,'aggressive',true,...
    'transitionmatrix',false,'stationarydistribution',[],...
    'teleportation',0.1); %set defaults
options.set(varargin); %set given options


W=sparse(W);
N=options.nodes;
aggressive=options.aggressive;

if ~options.isset('stationarydistribution')
    % check W is connected
    [WC,original_node_index]=LCC(W);
    
    if length(WC)~=length(W)
        warning('considering only largest connnected component');
        W=WC;
        N=min(N,length(WC));
    end
    clear('WC');
else
    % remove nodes with 0 mass in stationary distribution
    original_node_index=find(options.stationarydistribution);
    if length(original_node_index)~=length(W)
        warning('considering only nodes with non-zero mass in stationary distribution');
        W=W(original_node_index,original_node_index);
        N=min(N,length(original_node_index));
        options.stationarydistribution=options.stationarydistribution(original_node_index);
    end
end

if ~options.transitionmatrix
    if sum(diag(W))
        warning('removed self-edges')
        W=W-diag(diag(W));
    end
    k=sum(W,1);
    vol=sum(k);
    [row,col,val]=find(W);
    P=sparse(row,col,val./k(col)');
else
    P=W;
    vol=length(P);
end

if ~options.isset('stationarydistribution')
    if ~isequal(W,W')
        d=page_rank(P,options.teleportation,sum(W,2));
        [row,col,val]=find(P);
        d=d.*vol;
        W=sparse(row,col,val.*d(col),size(P,1),size(P,2));
    else
        d=sum(W,2);
    end
else
    d=options.stationarydistribution(:);
    [row,col,val]=find(P);
    d=d.*vol;
    W=sparse(row,col,val.*d(col),size(P,1),size(P,2));
end


if options.isset('local')
    local=options.local;
    %reindex to LCC
    if iscell(local)
        for i=1:length(local)
            local{i}=find(ismember(original_node_index,local{i}));
        end
        local=local(~cellfun(@isempty,local));
    else
        local=find(ismember(original_node_index,local));
    end
    if length(local)>N
        node_sampler=@() local(randsample(1:length(local),N,false));
    else
        node_sampler=@() local(randperm(length(local)));
    end
else
    node_sampler=@() randsample(1:length(W),N,false);
end

min_viscount=options.viscount;

% set defaults for truncation and alpha based on cut_function
switch cut_function
    
    case 'MOV'
        f_handle=@MOVcut;
        
        if isempty(options.alpha)
            e=laplace_eig(W);
            a=linspace(0.7,1/(1-e),20);
            alpha=(a-1)./a-10^-10;
        else
            alpha=options.alpha;
        end
        if isempty(options.truncation)
            truncation=inf;
        else
            truncation=options.truncation;
        end
        
    case 'ACL'
        f_handle=@ACLcut;
        
        if isempty(options.alpha)
            alpha=.001;
        else
            alpha=options.alpha;
        end
        if isempty(options.truncation)
            nmax=full(max(sum(W)));
            m=full(sum(sum(W)));
            truncation=logspace(-log10(nmax),-log10(m),20);
        else
            truncation=options.truncation;
        end
        
    case 'EGO'
        f_handle=@EGOcut;
        alpha=1;
        truncation=1;
        
    otherwise
        error('NCP:cut_function','unknown cut function');
end

%preallocate
conductance_con=ones(length(W)-1,1)*inf;
if nargout>1
    communities_con=cell(length(W)-1,1);
end
if nargout>2
    conductance_dis=ones(length(W)-1,1)*inf;
end
if nargout>3
    communities_dis=cell(length(W)-1,1);
end
if nargout>4
    assoc_mat=zeros(length(W));
end

for j=1:length(alpha)
    for k=1:length(truncation)
        visitcount=zeros(length(W),1);
        i=1;
        disp(['alpha =',num2str(alpha(j)),', trunc = ',num2str(truncation(k))]);
        %random reordering of nodes to avoid always sampling same subset
        %due to visitcount reached
        nodes=node_sampler();
        while min(visitcount)<min_viscount&&i<=length(nodes)
            
            [supp,cond,flag,connected]=f_handle(W,d,nodes(i),alpha(j),truncation(k));
            
            if flag
                warning('maximum number of iterations reached while computing pagerank, skipping to next parameter value')
                break
            end
            
            if ~isempty(supp)
                [~,mk]=min(cond);
                visitcount(supp(1:mk))=visitcount(supp(1:mk))+1;
            end
            
            n_supp=length(supp);
            update=find((cond(:)<conductance_con(1:n_supp))&connected(:));
            conductance_con(update)=cond(update);
            
            if nargout>1
                for l=update(:)'
                    communities_con{l}=original_node_index(supp(1:l));
                end
            end
            
            if nargout>2
                update_dis=find(cond(:)<conductance_dis(1:n_supp));
                conductance_dis(update_dis)=cond(update_dis);
                if nargout>3
                    for l=update(:)'
                        communities_dis{l}=original_node_index(supp(1:l));
                    end
                end
            end
            
            if nargout>4
                if ~isempty(supp)
                    assoc_mat(supp(1:mk),supp(1:mk))=assoc_mat(supp(1:mk),supp(1:mk))+1;
                end
            end
            
            if min(visitcount)>=min_viscount
                disp(['min_viscount reached with ',num2str(min(visitcount)),' visits after ',num2str(i),' iterations']);
            end
            
            if aggressive
                if iscell(nodes)
                    remove=find(cellfun(@(nn) all(visitcount(nn)>=min_viscount),nodes(i+1:end)));
                else
                    remove=find(visitcount(nodes(i+1:end))>=min_viscount);
                end
                nodes((i)+remove)=[];
            end
            
            i=i+1;
        end
        fprintf('sampled %u nodes\n',i-1);
    end
end

end
