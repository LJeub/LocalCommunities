function [conductance_con,communities_con,conductance_dis,communities_dis,assoc_mat]=multilayerNCP(W,cut_function,varargin)

options=OptionStruct('walktype','classical','layercoupling',1,'teleportation',0);
NCPoptions=OptionStruct('nodes',length(W)*length(W{1}),'local',[],'alpha',[],'truncation',[],...
    'viscount',10,'aggressive',true,'transitionmatrix',false,'stationarydistribution',[],'teleportation',[]);
ncpopts=options.setvalid(varargin);
NCPoptions.set(ncpopts);
if ~isstruct(NCPoptions)
    NCPoptions=NCPoptions(:);
    NCPoptions=struct(NCPoptions(1:2:end),NCPoptions(2:2:end),1);
end

switch options.walktype
    case 'classical'
        W=supra_adjacency(W,options.layercoupling);
        if options.teleportation>0
            kin=sum(W,2);
            kout=sum(W,1);
            W=W*diag(kout.^-1);
            p=page_rank(W,options.teleportation,kin);
            NCPoptions.stationarydistribution=p;
            NCPoptions.transitionmatrix=true;
        end
    case 'relaxed'
        P=relax_rate_walk(W);
        W=spblkdiag(W{:});
        kin=sum(W,2);
        W=P(options.layercoupling);
        p=page_rank(W,options.teleportation,kin);
        NCPoptions.stationarydistribution=p;
        NCPoptions.transistionmatrix=true;
end

switch nargout
    case {0,1}
        [conductance_con]=NCP(W,cut_function,NCPoptions);
    case 2
        [conductance_con,communities_con]=multilayerNCP(W,cut_function,NCPoptions);
    case 3
        [conductance_con,communities_con,conductance_dis]=multilayerNCP(W,cut_function,NCPoptions);
    case 4
        [conductance_con,communities_con,conductance_dis,communities_dis]=multilayerNCP(W,cut_function,NCPoptions);
    case 5
        [conductance_con,communities_con,conductance_dis,communities_dis,assoc_mat]=multilayerNCP(W,cut_function,NCPoptions);
end
        
end

