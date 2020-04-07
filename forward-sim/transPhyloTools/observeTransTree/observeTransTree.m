function [observedNodes, observedLineList] = observeTransTree(config)
% function to sample linelist struct based on metadata
% 
% INPUTS:
%   observeConfig: struct array with fields
%       files:
%           infectionLinelistFilename : filename of linelist data
%           [OPTIONAL] outputFilename : filename for to output list of
%              sampled nodes
%       samplingFrame:
%           type: 
%           levels: cellarray of cellarrays with levels for each factor
%       probabilityModel: struct defining probability sampling model,
%          possibly with nested covariates
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


%% Sample linelist

% top level eligibility
if strcmp(config.samplingFrame.probabilityModel.type,'weightedRandomSample')
    if isfield(config.samplingFrame.probabilityModel,'timeInterval')
        % prune linelist of ineligible times
        ineligibleTimeIntervalIdx = ([linelist.timeInfected] < config.samplingFrame.probabilityModel.timeInterval{1}) | ([linelist.timeInfected] > config.samplingFrame.probabilityModel.timeInterval{2});
        linelist = linelist(~ineligibleTimeIntervalIdx);
    end
    %prune off ineligible levels
    ix=ismember({linelist.(config.samplingFrame.type)},config.samplingFrame.levels);
    linelist=linelist(ix);
else
    error('only weightedRandomSample is implemented')
end

% second level elegibility
secondaryCategories=fieldnames(config.samplingFrame);
secondaryCategories=secondaryCategories(~ismember(secondaryCategories,{'type','levels','probabilityModel'}));
if ~isempty(secondaryCategories)
    for n=1:length(secondaryCategories)
        % prune linelist of ineligible times
        if isfield(config.samplingFrame.(secondaryCategories{n}).probabilityModel,'timeInterval')
            inCategoryIdx = ismember({linelist.(config.samplingFrame.type)},secondaryCategories(n));
            ineligibleTimeIntervalIdx = ([linelist.timeInfected] < config.samplingFrame.(secondaryCategories{n}).probabilityModel.timeInterval{1}) | ([linelist.timeInfected] > config.samplingFrame.(secondaryCategories{n}).probabilityModel.timeInterval{2});
            linelist = linelist(~(inCategoryIdx & ineligibleTimeIntervalIdx));
        end
        %prune off ineligible levels
        if config.samplingFrame.(secondaryCategories{n}).probabilityModel.isExclusive == true
            ix=ismember({linelist.(config.samplingFrame.(secondaryCategories{n}).type)},config.samplingFrame.(secondaryCategories{n}).levels);
            linelist=linelist(ix);
        end
    end
end

% top level probability
[~,loc]=ismember({linelist.(config.samplingFrame.type)},config.samplingFrame.levels);

    % find probability of level in sample
    weights=histcounts(loc).*config.samplingFrame.probabilityModel.probability;
    
    % sort to secondaryCategroies order
    if ~isempty(secondaryCategories)
        [ix,loc]=ismember(secondaryCategories,config.samplingFrame.levels);
        weights=weights(loc(ix));
    end
    weights=weights/sum(weights);

% list of lists of categories and probablities
topLevel={};
secondLevel ={};
samplingProbability=[];
if ~isempty(secondaryCategories)
    for n=1:length(secondaryCategories)
        % assembly total sampling probabiltiy for each class of subject
        for m=1:length(config.samplingFrame.(secondaryCategories{n}).levels)
            topLevel{n,m}=[{config.samplingFrame.type},secondaryCategories(n)];
            secondLevel{n,m}=[config.samplingFrame.(secondaryCategories{n}).type,config.samplingFrame.(secondaryCategories{n}).levels(m)];
            samplingProbability(n,m) = weights(n)*config.samplingFrame.(secondaryCategories{n}).probabilityModel.probability{m};
        end
    end
else
    for k=1:length(config.samplingFrame.levels)
        topLevel{k}=[{config.samplingFrame.type},config.samplingFrame.levels(k)];
    end
    samplingProbability=weights;
end
% flatten
topLevel={topLevel{:}};
secondLevel={secondLevel{:}};
samplingProbability=cumsum([samplingProbability(:)]);
samplingProbability = samplingProbability/max(samplingProbability);
keepIdx=~cellfun(@isempty,topLevel);
topLevel=topLevel(keepIdx);
if ~isempty(secondLevel)
    secondLevel=secondLevel(keepIdx);
end
samplingProbability=samplingProbability(keepIdx);

% random sample without replacement
observedNodes=cell(config.dataRequested.numSamples,1);
for k=1:config.dataRequested.numSamples
    
    valid=false;
    while ~valid
        samplingBin = 1+sum(samplingProbability<rand);
        eligibleIdx = ismember({linelist.(topLevel{samplingBin}{1})},topLevel{samplingBin}(2));
        if ~isempty(secondLevel)
                eligibleIdx = eligibleIdx & ismember({linelist.(secondLevel{samplingBin}{1})},secondLevel{samplingBin}(2));
        end
        valid = sum(eligibleIdx)>0;
    end
    
    node = randsample({linelist(eligibleIdx).id},1);
    count=0;
    while count<20 && k>1 && ismember(node,observedNodes(1:k-1)) || ~isempty(regexp(node{1}, 'ROOT', 'once')) 
        node = randsample({linelist(eligibleIdx).id},1);
        count=count+1;
        if count==20 % redraw
            warning(['could not find eligble sample in ',topLevel{samplingBin}{2},'-',secondLevel{samplingBin}{2}]);
             valid=false;
             while ~valid
                samplingBin = 1+sum(samplingProbability<rand);
                eligibleIdx = ismember({linelist.(topLevel{samplingBin}{1})},topLevel{samplingBin}(2));
                if ~isempty(secondLevel)
                        eligibleIdx = eligibleIdx & ismember({linelist.(secondLevel{samplingBin}{1})},secondLevel{samplingBin}(2));
                end
                valid = sum(eligibleIdx)>0;
            end
            count=0;
        end
    end
    observedNodes(k)=node;
    
    if mod(100*k/length(observedNodes),1)<100*(1/length(observedNodes))
       display(['percent complete: ',num2str(k/length(observedNodes)*100)]);
    end
end


%% observed linelist
observedIdx=ismember({linelist.id},observedNodes);
observedLineList=linelist(observedIdx);
if isfield(config.dataRequested,'fields')
    notRequested = setdiff(fieldnames(observedLineList),config.dataRequested.fields);
    for k=1:length(notRequested)
        observedLineList=rmfield(observedLineList,notRequested{k});
    end
end

%% write out
if isfield(config.files,'observedLineListFilename') 
    writetable(struct2table(observedLineList), config.files.observedLineListFilename)
end

end