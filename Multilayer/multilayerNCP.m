function [conductance_con,communities_con,conductance_dis,communities_dis,assoc_mat]=multilayerNCP(A,cut_function,varargin)
% Convenience wrapper around NCP for multiplex networks
%
% 
%
% Input:
%           A: cell array of adjacency matrices for each layer of the
%           network (all layers need to be the same size)
%
%           cut_function: 'ACL','MOV','EGO' to select algorithm to identify
%           communities
%
%           options:
%                       walktype: 'classical' or 'relaxed'
%                           [default: 'classical']
%
%                       layercoupling: strength of interlayer edges for
%                           'classical' walk or relax rate for 'relaxed' walk
%                           [default: 1]
%
%                       teleportation: teleportation rate (useful for
%                           directed networks) for unrecorded link
%                           teleportation
%                           [default: 0]
%
%                       physicalnodes: when set to true, sample NCP
%                           using physical nodes
%                           [default: false]
%
%                       + all options for NCP
%
%                       Note about 'local' option for NCP:
%                       If option 'physicalnodes' is set to false,
%                       nodes can be specified either using the state-node
%                       id (single number from 1:#nodes*#layers) or as a
%                       pair of node id and layer id. If providing two state
%                       indeces make sure they are provided as a column
%                       vector (otherwise it's treated as a node-layer
%                       pair). Node-layer pairs need to be provided as a
%                       nx2 matrix, where each row is a node-layer pair.
%
%                       If option 'physicalnodes' is set to true, 'local'
%                       should contain indeces of physical nodes.
%
% Output:
%           conductance_con: vector of minimum conductance values
%               considering only connected communities, e.g.
%               conductance_con(10) gives the minimum conductance found for
%               communities with 10 state nodes that are connected.
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

% Version: 1.02
% Date: Fri 18 Jul 2014 13:16:55 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

options=OptionStruct('walktype','classical','layercoupling',1,'teleportation',0,'physicalnodes',false);
NCPoptions=OptionStruct('nodes',length(A)*length(A{1}),'local',[],'alpha',[],'truncation',[],...
    'viscount',10,'aggressive',true,'transitionmatrix',false,'stationarydistribution',[],'teleportation',[]);
ncpopts=options.setvalid(varargin);
NCPoptions.set(ncpopts);

N=length(A{1});
% deal with physicalnodes option by setting the local option for NCP
if options.physicalnodes
    p_nodes=cell(N,1);
    for i=1:N
        p_nodes{i}=[repmat(i,length(A),1),(1:length(A))'];
    end
    if NCPoptions.isset('local')
        if iscell(NCPoptions.local)
            local=cell(length(NCPoptions.local),1);
            for i=1:length(NCPoptions.local)
                local{i}=vertcat(p_nodes{NCPoptions.local{i}});
            end
            NCPoptions.local=local;
        else
            NCPoptions.local=p_nodes(NCPoptions.local);
        end
    else
        NCPoptions.local=p_nodes;
    end
end

% convert 'local' option given as nodelayer indeces to state indeces
if NCPoptions.isset('local')
    NCPoptions.local=nodelayer2state(N,NCPoptions.local);
end

% set up different walk types
switch options.walktype
    case 'classical'
        A=supra_adjacency(A,options.layercoupling);
        kin=sum(A,2);
        kout=sum(A,1);
        A=A*diag(kout.^-1);
        p=page_rank(A,options.teleportation,kin);
        if any(kin==0)
            if max(p(kin==0))<10^-10
                p(kin==0)=0; % ensure these get removed later
            else
                error('excessive numerical error in page_rank calculation')
            end
        end
        NCPoptions.stationarydistribution=p;
        NCPoptions.transitionmatrix=true;
        
    case 'relaxed'
        P=relax_rate_walk(A);
        A=spblkdiag(A{:});
        kin=sum(A,2);
        A=P(options.layercoupling);
        p=page_rank(A,options.teleportation,kin);
        if any(kin==0)
            if max(p(kin==0))<10^-10
                p(kin==0)=0; % ensure these get removed later
            else
                error('excessive numerical error in page_rank calculation')
            end
        end
        NCPoptions.stationarydistribution=p;
        NCPoptions.transitionmatrix=true;
end

% Call NCP with appropriate number of outputs for efficiency
switch nargout
    case {0,1}
        [conductance_con]=NCP(A,cut_function,NCPoptions);
    case 2
        [conductance_con,communities_con]=NCP(A,cut_function,NCPoptions);
        communities_con=state2nodelayer(N,communities_con);
    case 3
        [conductance_con,communities_con,conductance_dis]=NCP(A,cut_function,NCPoptions);
        communities_con=state2nodelayer(N,communities_con);
    case 4
        [conductance_con,communities_con,conductance_dis,communities_dis]=NCP(A,cut_function,NCPoptions);
        communities_con=state2nodelayer(N,communities_con);
        communities_dis=state2nodelayer(N,communities_dis);
    case 5
        [conductance_con,communities_con,conductance_dis,communities_dis,assoc_mat]=NCP(A,cut_function,NCPoptions);
        communities_con=state2nodelayer(N,communities_con);
        communities_dis=state2nodelayer(N,communities_dis);
end

end



