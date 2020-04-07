function linelist = addImportationsField(config)
% function parses transmission tree to find independent roots, and labels
% them as importations
%
% INPUTS:
%   spatialConfig: struct array with fields
%       files:
%           infectionLinelistFilename : filename of linelist data
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

transmissionTree = transTreeFromLineList(linelist);

% find connected components
bin = conncomp(transmissionTree,'Type','weak');
rootNodes=cell(1,max(bin));
for k=1:max(bin)
    idx = bin == k;
    SG = subgraph(transmissionTree, idx);
    node = SG.Nodes.Name{1};
    predNode = predecessors(SG,node);
    while ~isempty(predNode) 
        node = predNode{1};
        predNode = predecessors(SG,node);
    end
    rootNodes{k}=node;
end
roots=rootNodes(bin);

% sort to linelist order
[ix,loc]=ismember({linelist.id},transmissionTree.Nodes.Name);

% populate
[linelist.importation]=roots{loc(ix)};


%% write out
if ~isfield(config.files,'outputFilename') % overwrite
    writetable(struct2table(linelist), config.files.infectionLinelistFilename)
else
    writetable(struct2table(linelist), config.files.outputFilename)
end

end