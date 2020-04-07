function linelist = addFactorAttributes(config)
% function to add factors to linelist struct
% factors like age, sex, vaccination status, hasFever, etc
% 
% INPUTS:
%   factorConfig: struct array with fields
%       files:
%           infectionLinelistFilename : filename of linelist data
%           [OPTIONAL] outputFilename : filename for output if you don't
%              want to overwrite
%       factors: cellarray of factor names
%       levels: cellarray of cellarrays with levels for each factor
%       probability: cellarray of vectors with multinomial probality of
%           assigning each factor
%       conditionalFactors: nested struct of factors that depend on top
%       level factors
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

%% add factors
% get primary fields
topNames=fieldnames(config);
topNames=topNames(~ismember(topNames,{'factors','files'}));

for k=1:length(topNames)
    secondNames{k}=fieldnames(config.(topNames{k}));
    secondNames{k}=secondNames{k}(~ismember(secondNames{k},{'levels','probability'}));
end
if ~all(ismember(config.factors,[topNames;secondNames{:}]))
    error('conditional factors can only be one level deep for now because Mike is bad at recursion.')
end

for k=1:length(topNames)
    cumProb=cumsum(cell2mat(config.(topNames{k}).probability))/sum(cell2mat(config.(topNames{k}).probability));
    factorData = config.(topNames{k}).levels(1+sum(cumProb<rand(length(linelist),1),2));
    [linelist.(topNames{k})]=factorData{:};
    for n=1:length(secondNames{k})
        factorData=cell(size(factorData));
        for m=1:length(config.(topNames{k}).levels)
            tmp=config.(topNames{k}).(secondNames{k}{n}).(config.(topNames{k}).levels{m});
            cumProb=cumsum(cell2mat(tmp.probability))/sum(cell2mat(tmp.probability));
            idx=ismember([linelist.(topNames{k})],config.(topNames{k}).levels{m});
            factorData(idx) = tmp.levels(1+sum(cumProb<rand(sum(idx),1),2));
        end
        [linelist.(secondNames{k}{n})]=factorData{:};
    end
end



%% write out
if ~isfield(config.files,'outputFilename') % overwrite
    writetable(struct2table(linelist), config.files.infectionLinelistFilename)
else
    writetable(struct2table(linelist), config.files.outputFilename)
end

end