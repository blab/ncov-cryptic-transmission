function totalInfections = treeClusterSim(config)
% protype of treeClusterSim

%% intialize 
if isstruct(config)
    
elseif ischar(config)
    config=yaml.ReadYaml(config);
else 
    error('config must be struct or filename');
end

if ~exist(['output/',config.simulation.name],'dir')
    mkdir(['output/',config.simulation.name]);
end

% simulation
    startTime=config.simulation.startTime;
    dt=config.simulation.timeStep;
    endTime=startTime+config.simulation.totalTime;
    
% rng seed
%     distribution(config, 'setSeed')


% internals
    global clusterIdList
    clusterIdList={};

    global infectionIdList;
    infectionIdList=repmat({''},1,1e5);
    
    global totalInfections;
    totalInfections=0;

% collect all clusters and infections to write out
    global clusterLineList;
    clusterLineList=cell(0,8);
    
    global infectionLineList;
    infectionLineList=cell(1e5,6);                    


    %% run

t=startTime+dt;
while t<endTime 
    
    %% initial condition 
    nClust=0;
    switch config.initialConditions.type
        case 'fixed'
            if t==(startTime+dt)
                switch config.clusterModel.type
                    case 'stochasticSEIR'
                        nClust=distribution(config.initialConditions.clusters);
                        for k=1:nClust
                            clusterObj(k)=cluster.stochasticSEIR(config,startTime,['ROOT_',num2str(k)]);
                            clusterObj(k).N=distribution(config.initialConditions.clusters.params.N);

                            nLatents = distribution(config.initialConditions.clusters.params.latents);
                            infectedById = cell(1,nLatents);
                            for n=1:nLatents
                                infectedById{n}=[clusterObj(k).id,'_ROOT'];
                            end
                            clusterObj(k).createNewInfections(nLatents,config,startTime,infectedById);

                        end
                    otherwise
                        error('unknown clusterModel type')
                end
            end
        case 'multiple'
            start_times = cell2mat(config.initialConditions.start_times)+dt;
            switch config.clusterModel.type
                case 'stochasticSEIR'
                    for k=1:length(start_times)
                        if start_times(k)>=t && start_times(k)<t+dt
                            if ~exist('clusterObj')
                                clusterObj=cluster.stochasticSEIR(config,t-dt,['ROOT_',num2str(1)]);
                                clusterObj.N=distribution(config.initialConditions.clusters.(['set_',num2str(k)]).params.N);
                            else
                                clusterObj(end+1)=cluster.stochasticSEIR(config,t-dt,['ROOT_',num2str(k)]);
                                clusterObj(end).N=distribution(config.initialConditions.clusters.(['set_',num2str(k)]).params.N);
                            end

                            nLatents = distribution(config.initialConditions.clusters.(['set_',num2str(k)]).params.latents);
                            infectedById = cell(1,nLatents);
                            for n=1:nLatents
                                infectedById{n}=[clusterObj(k).id,'_ROOT'];
                            end
                            clusterObj(k).createNewInfections(nLatents,config,startTime,infectedById);

                            nClust=nClust +1;
                        end
                    end
                otherwise
                    error('unknown clusterModel type')
            end

        otherwise
            error('unknown initial conditions type')
    end
    

    
    % clean up from last round
    killList=false(length(clusterObj),1);
    for k=1:length(clusterObj)
        clusterObj(k).removeClearedInfections; 
        if clusterObj(k).currentInfections==0
            killList(k)=true;
            cluster.stochasticSEIR.cacheLineList(clusterObj(k),t);
        end
    end
    clusterObj=clusterObj(~killList);
    
    % transmit within each cluster
    switch config.clusterModel.type
        case 'stochasticSEIR'

            for k=1:length(clusterObj)
                clusterObj(k).transmit(config,t);
            end     
    end
    
    % spawn new clusters from existing
    switch config.branchingModel.type
        case 'perInfection'
            lambdaNew = [clusterObj.newClusterBeta].*[clusterObj.currentInfections]*config.clusterModel.params.infectiousDuration.params.mu;
            for k=1:length(lambdaNew)
                newC=poissrnd(lambdaNew(k));
                for n=1:newC
                    infectedId={clusterObj(k).infections.id};
                    infectedId=infectedId(randi(length(infectedId),1));

                    clusterObj(end+1)=cluster.stochasticSEIR(config,t,clusterObj(k).id);
                    clusterObj(end).createNewInfections(1,config,t,infectedId);
                end
            end
        case 'perClusterPerDay'
            lambdaNew = [clusterObj.newClusterBeta].*ones(length(clusterObj)/1)*dt;
            for k=1:length(lambdaNew)
                newC=poissrnd(lambdaNew(k));
                for n=1:newC
                    infectedId={clusterObj(k).infections.id};
                    infectedId=infectedId(randi(length(infectedId),1));

                    clusterObj(end+1)=cluster.stochasticSEIR(config,t,clusterObj(k).id);
                    clusterObj(end).createNewInfections(1,config,t,infectedId);
                end
            end    
    end
    
    % spawn new clusters from importation
    if isfield(config,'importationModel')
        switch config.importationModel.type
            case 'denovoPerTimestep'
                beta = distribution(config.importationModel.params.newClusterBeta);
                
                newC = poissrnd(beta*dt);

                for k=1:newC

                    nClust = nClust+1;
                    sourceCluster=['ROOT_',num2str(nClust)];

                    clusterObj(end+1)=cluster.stochasticSEIR(config,t,sourceCluster);

                    clusterObj(end).N=distribution(config.importationModel.params.N);

                    nInfections = distribution(config.importationModel.params.infections);
                    infectedById = cell(1,nInfections);
                    for n=1:nInfections
                        infectedById{n}=[clusterObj(end).id,'_ROOT'];
                    end
                    clusterObj(end).createNewInfections(nInfections,config,t,infectedById);

                end
        end
    end
    
    % quit if 0 infections and no importations
    if ~isfield(config,'importationModel') && isempty(clusterObj)
        t=endTime;
    end
    
    t=t+dt;
    if mod(t,(endTime-startTime)/10)<dt
        t
        totalInfections
    end
    
end
t
totalInfections

% clean up from last round
    for k=1:length(clusterObj)
        clusterObj(k).removeClearedInfections; 
        cluster.stochasticSEIR.cacheLineList(clusterObj(k),t);
    end
    

%% write out line list data
cluster.stochasticSEIR.writeCachedLineList(config);
infection.SEIR.writeCachedLineList(config);

end