classdef stochasticSEIR < handle
   properties %(SetAccess = immutable) 
       id = '';
       infections;
       totalInfectedIds=[];
       currentInfections=nan;
       newClusterBeta = nan;
       seededById = nan;
       timeCreated = nan;
       N = nan;
   end
   methods 
      function obj = stochasticSEIR(config,timeCreated,seededById)
          if nargin > 0
             obj = obj.uid();
             obj.currentInfections=0;
             obj.newClusterBeta = distribution(config.branchingModel.params.branchingBeta);
             obj.N = ceil(distribution(config.clusterModel.params.N));
          end
          if nargin>1
              obj.timeCreated=timeCreated;
          end
          if nargin>2
              obj.seededById = seededById;
          end
      end
      
      function obj = uid(obj)

          global clusterIdList
 
          str='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
          maxSize=7;
                     
          id=str(randi(maxSize,1,maxSize));
          if any(strcmp(id,clusterIdList))
              id=[id,str(randi(maxSize,1,1))]; % make uid unique 
          end
          
          clusterIdList{end+1}=id;
          obj.id=id;
      end
      
      function obj = createNewInfections(obj,newInfections,config,t,infectedById)
          for k=1:newInfections
             if isempty(obj.infections)
                obj.infections=infection.SEIR(config,t,infectedById{k},obj.id);
             else
                obj.infections(end+1)=infection.SEIR(config,t,infectedById{k},obj.id);
             end
          end
          obj.totalInfectedIds=union(obj.totalInfectedIds,{obj.infections.id});
          obj.currentInfections = length([obj.infections]);
      end

      
      function obj = removeClearedInfections(obj)
          killList=false(length(obj.infections),1);
          for k=1:length(killList)
              if obj.infections(k).infectedDuration >= obj.infections(k).latentDuration + obj.infections(k).infectiousDuration
                  killList(k)=true;
              end
          end
          obj.infections=obj.infections(~killList);
          obj.currentInfections = length(obj.infections);
      end
      
      function obj = transmit(obj,config,t)
          if ~isempty(obj.infections)
            obj.infections.updateInfection(config);

            fS=1-length(obj.totalInfectedIds)/obj.N;
            amp = 1 + config.clusterModel.params.seasonality.params.amplitude * cos(2*pi/config.clusterModel.params.seasonality.params.period*(t-config.clusterModel.params.seasonality.params.phase));
            pTrans=1-exp(-[obj.infections.beta].*double([obj.infections.sheddingDuration]>0).*fS * amp * config.simulation.timeStep);
            sourceUID={obj.infections.id};

            newI = pTrans>=rand(size(pTrans));

            if any(newI)
                obj.createNewInfections(sum(newI),config,t,sourceUID(newI));
            end
            obj.currentInfections = length(obj.infections);
          end
      end
   end
   methods(Static)
      function cacheLineList(obj,t)
          % has to be done on cluster "death" because I want total infection history data
          global clusterLineList
          
          clusterLineList=[clusterLineList; [{obj.id},{obj.timeCreated},{t},{obj.N}, ...
              {length(obj.totalInfectedIds)},{obj.currentInfections},{obj.newClusterBeta},{obj.seededById}]];
          
      end
      function writeCachedLineList(config)
          global clusterLineList

          fid = fopen(config.files.clusterLinelistFilename,'w+');
          fprintf(fid,'id,timeCreated,timeKilled,N,totalInfections,currentInfections,newClusterBeta,seededById\n');
          for k=1:size(clusterLineList,1)
              fprintf(fid,'%s,%f,%f,%f,%f,%f,%f,%s\n' ...     
                  ,clusterLineList{k,1},clusterLineList{k,2}, ...
                  clusterLineList{k,3},clusterLineList{k,4},clusterLineList{k,5}, ...
                  clusterLineList{k,6},clusterLineList{k,7},clusterLineList{k,8});
          end
          fclose(fid);
      end
   end
end