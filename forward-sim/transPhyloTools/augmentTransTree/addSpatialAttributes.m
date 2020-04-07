function linelist = addSpatialAttributes(config)
% function to add spatial location to linelist
% 
% INPUTS:
%   spatialConfig: struct array with fields
%       files:
%           infectionLinelistFilename : filename of linelist data
%       location: cellarray of names describing locations (ie GEOID)
%       connectivityGraph: cellarray of cellarrays with row of indexes for
%           which locations are linked (adjacency matrix of GMRF precision)
%   
%   OPTIONAL:
%       hoppingModel: config that defines model of hopping on (and off) the
%       connectivityGraph
% 
%       OR
%
%       [NOT YET IMPLEMENTED] connectivityWeights: cellarray of cellarrays with weights for linked
%           locations defined by connectivityGraph
%
%       -------------
%       aggregation:
%           level1: array of names of higher level of location.name
%           level2: same for another aggregation
%           ...
%   
%  OUTPUT:
%   linelist returned with new columns defined by factors
%

%% initialize 
if isstruct(config)
    
elseif ischar(config)
    config=yaml.ReadYaml(config);
else 
    error('config must be struct or filename');
end

linelist=xls2struct(config.files.infectionLinelistFilename,'structArray');

transTree = transTreeFromLineList(linelist);

% clean up config
config.location.rowId=cell2mat(config.location.rowId);
config.location.name=[config.location.name{:}];
for k=1:length(config.connectivityGraph)
    config.connectivityGraph{k}=cell2mat(config.connectivityGraph{k});
end

edgeList=transTree.Edges{:,1};

%% add random root locations
rootNodes=transTree.Nodes.Name(~cellfun(@isempty,regexp(transTree.Nodes.Name,'ROOT')));
rootLoc = randsample(config.location.rowId,length(rootNodes));

%% run the spatial model and add the attributes to the linelist

lookBack=100;
newLocs=cell(height(transTree.Edges),1);
for k=1:length(newLocs)
    
    ancNode=edgeList(k,1);
    
    % if new importation, start from root
    if ~cellfun(@isempty,regexpi(ancNode,'ROOT'))
        locIdx= rootLoc(ismember(rootNodes,ancNode));
        newLocs{k}=config.location.name(locIdx);
    else
        ancLocIdx = max(0,k-lookBack)+find(ismember(edgeList(max(1,k-lookBack+1):k,2),ancNode));
        if isempty(ancLocIdx)
            ancLocIdx = find(ismember(edgeList(1:k,2),ancNode));
            lookBack = k-ancLocIdx+1;
        end
        ancLoc=newLocs{ancLocIdx};
        ancIdx=find(ismember(config.location.name,ancLoc));
        connectedLocs=config.connectivityGraph{ancIdx};
        
        % find transition probabilities
        if isfield(config,'hoppingModel')
            hopProb=nan(size(connectedLocs));
        
            % self
            hopProb(ancIdx==connectedLocs)=config.hoppingModel.selfProbability;
            
            % neighbors
            if strcmp(config.hoppingModel.neighborProbabilityModel,'oneOverN')
                neighborHopProb=(1-config.hoppingModel.selfProbability-config.hoppingModel.hopAnywhereProbability)/(length(hopProb)-1);
                hopProb(ancIdx~=connectedLocs)=neighborHopProb;
            else
                error('only oneOverN implemented currently.')
            end
            
            % uniformRandom
            if config.hoppingModel.hopAnywhereProbability>0
                hopProb(end+1)=config.hoppingModel.hopAnywhereProbability;
            end
            hopProb=hopProb/sum(hopProb);
            
        else
            error('only hoppingModel is currently implemented.')
        end
        
        % select new location
        locIdx=1+sum(cumsum(hopProb)<=rand);
        if config.hoppingModel.hopAnywhereProbability>0 && locIdx==length(hopProb)
            randLoc=randsample(config.location.name,1);
            newLocs{k}=randLoc;
        else
            newLocs{k}=config.location.name(connectedLocs(locIdx));
        end

    end
    if mod(100*k/length(newLocs),1)<100*(1/length(newLocs))
       display(['percent complete: ',num2str(k/length(newLocs)*100)]);
    end
end

%% add to linelist 
[~,loc]=ismember({linelist.id}, edgeList(:,2));
[linelist.(config.location.type)]=deal(newLocs{loc});


%% add aggregated levels

% for fields in aggregation
% look up indices of each level
% map level and output

if isfield(config,'aggregation')
    rowIdx=cellfun(@(x) find(ismember(config.location.name,x)),{linelist.(config.location.type)});
    names=fieldnames(config.aggregation);
    for k=1:length(names)
        [linelist.(names{k})]=config.aggregation.(names{k}){rowIdx};
    end    
end

%% write out
if ~isfield(config.files,'outputFilename') % overwrite
    writetable(struct2table(linelist), config.files.infectionLinelistFilename)
else
    writetable(struct2table(linelist), config.files.outputFilename)
end


end