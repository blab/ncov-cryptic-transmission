classdef SEIR < handle
   properties %(SetAccess = immutable) 
      id = '';
      infectedDuration = 0;
      latentDuration = nan;
      infectiousDuration = nan;
      sheddingDuration = 0;
      infectedById = '';
      timeInfected = nan;
      beta = nan;  
   end
   methods 
      function obj = SEIR(config,timeInfected,infectedById,clusterId)
          if nargin > 0
             if isnumeric(config.clusterModel.params.infectiousDuration)
                obj.infectiousDuration = config.clusterModel.params.infectiousDuration;
             elseif isstruct(config.clusterModel.params.infectiousDuration)
                 obj.infectiousDuration = distribution(config.clusterModel.params.infectiousDuration);
             else
                error('Value must be numeric or param struct')
             end
             if isnumeric(config.clusterModel.params.latentDuration)
                obj.latentDuration = config.clusterModel.params.latentDuration;
             elseif isstruct(config.clusterModel.params.latentDuration)
                 obj.latentDuration = distribution(config.clusterModel.params.latentDuration);
             else
                error('Value must be numeric or param struct')
             end
             if isnumeric(config.clusterModel.params.beta)
                obj.beta = config.clusterModel.params.beta;
             elseif isstruct(config.clusterModel.params.beta)
                 obj.beta = distribution(config.clusterModel.params.beta);
             else
                error('Value must be numeric or param struct')
             end
             if nargin>1
                 obj.timeInfected=timeInfected;
             end
             if nargin>2
                 obj.infectedById = infectedById;
             end
             
             obj = obj.uid(clusterId);
          end
      end
      function obj = uid(obj,clusterId)

          global infectionIdList
          global totalInfections
          global infectionLineList
 
          if totalInfections==length(infectionIdList)
              infectionIdList=[infectionIdList,repmat({''},1,1e5)];
              infectionLineList=[infectionLineList;cell(1e5,6)];
          end
          totalInfections=totalInfections+1;
          
          obj.id=[clusterId,'_',num2str(totalInfections)];
          
          infectionIdList{totalInfections}=obj.id;
          
          infectionLineList(totalInfections,:)=[{obj.id},{obj.timeInfected},{obj.beta}, ...
              {obj.latentDuration},{obj.infectiousDuration},{obj.infectedById}];
      end
      function obj = updateInfection(obj,config)
          for k=1:length(obj)
            obj(k).infectedDuration=obj(k).infectedDuration + config.simulation.timeStep;
            if obj(k).infectedDuration >= obj(k).latentDuration
                obj(k).sheddingDuration = obj(k).sheddingDuration + config.simulation.timeStep;
            end
          end
      end
   end
   methods(Static)
      function writeCachedLineList(config)
          global infectionLineList
          global totalInfections
          
          fid = fopen(config.files.infectionLinelistFilename,'w+');
          fprintf(fid,'id,timeInfected,beta,latentDuration,infectiousDuration,infectedById\n');
          for k=1:totalInfections
              fprintf(fid,'%s,%f,%f,%f,%f,%s\n' ...     
                  ,infectionLineList{k,1},infectionLineList{k,2}, ...
                  infectionLineList{k,3},infectionLineList{k,4},infectionLineList{k,5},infectionLineList{k,6});
          end
          fclose(fid);
      end
   end
end