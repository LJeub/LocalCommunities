function [conductance_con,communities_con,conductance_dis,communities_dis,assoc_mat]=NCP(W,cut_function,varargin)
% [conductance_con,communities_con,conductance_dis,communities_dis,assoc_mat]=NCP(W,cut_function,options)
%
% NCP:  Compute the Network Community Profile for a network
%
% This function (with the appropriate options listed below) can be used
% to compute approximate local and global NCPs for a network using any
% of the three methods in:
% Jeub, L. G. S., Balachandran, P., Porter, M. A., Mucha, P. J., 
% & Mahoney, M. W. (2014).
% Think Locally, Act Locally: The Detection of Small, Medium-Sized, and 
% Large Communities in Large Networks. arXiv:1403.3795 [cs.SI]
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
% For example, to compute an NCP for a network with adjacency matrix W
% using the ACLcut method with alpha=0.01, one would call the function as
%   conductance=NCP(W,'ACL','alpha',0.01);
%
% One can then plot the NCP using 
%   loglog(conductance)
%
% see also ACLcut MOVcut normalize_assoc_mat

% Version: 1.01
% Date: Tue 25 Mar 2014 17:22:32 GMT
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

% Parse Options
options=OptionStruct('nodes',length(W),'local',[],'alpha',[],'truncation',[],...
    'viscount',10,'transition_matrix',false); %set defaults
options.set(varargin); %set given options

% check W is symmetric, connected, no self-edges (otherwise fix)
W=sparse(W);

[WC,ind]=LCC(W);
N=options.nodes;
if length(WC)~=length(W)
    warning('considering only largest connnected component');
    W=WC;
    N=min(N,length(WC));
end
clear('WC');

if ~options.transition_matrix
    if sum(diag(W))
        warning('removed self-edges')
        W=W-diag(diag(W));
    end
    P=W./repmat(sum(W,1),size(W,1),1);
else
    P=W;
end

if ~isequal(W,W')
    [d,eigval]=eigs(P,1);
    if norm(eigval-1)>10^-10
        error('largest eigenvalue of transition matrix not equal to 1')
    end
    d=d/sum(d);
else
    d=sum(W,2);
    d=d/sum(d);
end


%renormalized network
%if options.transition_matrix
    W=P.*repmat(d(:)',size(W,1),1);
%end



if options.isset('local')
    local=options.local;
    if length(local)>N
        node_sampler=@() randsample(local,N,false);
    else
        node_sampler=@() local(randperm(length(local)));
    end
else
    node_sampler=@() randsample(1:length(W),N,false);
end

min_viscount=options.viscount;


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
            nmax=max(sum(W));
            m=sum(W(:));
            epsstep=(log2(m)-log2(nmax))/20;
            truncation=1./2.^(log2(nmax):epsstep:log2(m));
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



%chunking for parallelisationn (try matlabpool to get size (fails if
%parallel toolbox is not installed))
try
    pool_size=matlabpool('size');
    if pool_size==0
        chunk=1;
        parallel=false;
    else
        parallel=true;
        chunk=pool_size*10;
    end
catch
    chunk=1;
    parallel=false;
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


argout=nargout;

for j=1:length(alpha)
    for k=1:length(truncation)
        visitcount=zeros(length(W),1);
        i=1;
        disp(['alpha =',num2str(alpha(j)),', trunc = ',num2str(truncation(k))]);
        %random reordering of nodes to avoid always sampling same subset
        %due to visitcount reached
        nodes=node_sampler();
        while min(visitcount)<min_viscount&&i<=length(nodes)
            if parallel
                loopsize=min(chunk-1,N-i);
                loopnodes=nodes(i:i+loopsize);
                alphal=alpha(j);
                truncl=truncation(k);
                try
                    parfor it=1:length(loopnodes)
                        [suppl{it},condl{it},flag,connectedl{it}]=f_handle(W,d,loopnodes(it),alphal,truncl);
                        if flag
                            error('not converged')
                        end
                    end
                catch err
                    
                    if ~isempty(strfind(err.message,'not converged'))
                        warning('maximum number of iterations reached while computing pagerank, skipping to next parameter value')
                        break
                    else
                        rethrow(err);
                    end
                end
                
                for it=1:length(suppl)
                    supp=suppl{it};
                    cond=condl{it};
                    connected=connectedl{it};
                    n_supp=length(supp);
                    if ~isempty(supp)
                        [~,mk]=min(cond);
                        visitcount(supp(1:mk))=visitcount(supp(1:mk))+1;
                    end
                    
                    update=find((cond(:)<conductance_con(1:n_supp))&connected(:));
                    conductance_con(update)=cond(update);
                    if argout>1
                        for l=update(:)'
                            communities_con{l}=supp(1:l);
                        end
                    end
                    
                    if argout>2
                        update_dis=find((cond(:)<conductance_dis(1:n_supp)));
                        conductance_dis(update_dis)=cond(update_dis);
                        if argout>3
                            for l=update(:)'
                                communities_dis{l}=supp(1:l);
                            end
                        end
                    end
                    
                    if argout>4
                        if ~isempty(supp)
                            assoc_mat(supp(1:mk),supp(1:mk))=assoc_mat(supp(1:mk),supp(1:mk))+1;
                        end
                    end
                end
                i=i+length(suppl);
            else
                [supp,cond,flag,connected]=f_handle(W,d,nodes(i),alpha(j),truncation(k));
                if flag
                    warning('maximum number of iterations reached while computing pagerank \n skipping to next parameter value')
                    break
                end
                if ~isempty(supp)
                    [~,mk]=min(cond);
                    visitcount(supp(1:mk))=visitcount(supp(1:mk))+1;
                end
                n_supp=length(supp);
                
                update=find((cond(:)<conductance_con(1:n_supp))&connected(:));
                conductance_con(update)=cond(update);
                if argout>1
                    for l=update(:)'
                        communities_con{l}=supp(1:l);
                    end
                end
                
                if argout>2
                    update_dis=find(cond(:)<conductance_dis(1:n_supp));
                    conductance_dis(update_dis)=cond(update_dis);
                    if argout>3
                        for l=update(:)'
                            communities_dis{l}=supp(1:l);
                        end
                    end
                end
                
                if argout>4
                    if ~isempty(supp)
                        assoc_mat(supp(1:mk),supp(1:mk))=assoc_mat(supp(1:mk),supp(1:mk))+1;
                    end
                end
                
                i=i+1;
            end
            if min(visitcount)>=min_viscount
                disp(['min_viscount reached with ',num2str(min(visitcount)),' visits after ',num2str(i),' iterations']);
            end
        end
    end
end
