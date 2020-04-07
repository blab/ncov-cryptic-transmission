% protype of treeClusterSim
% v 0.01

clear; clc;

% include
addpath(genpath('../../transPhyloTools'));

%% run model
config=yaml.ReadYaml('config.yaml');
treeClusterSim(config)
disp('done!')


plotTreeClusterSim(config);

%% construct downsampled phylogenetic tree



infection=xls2struct(config.files.infectionLinelistFilename,'structArray');

transmissionTree = transTreeFromLineList(infection);

% sample nodes
    rng(1000)
    sampleIdx=(1/100*ones(height(transmissionTree.Nodes),1))>rand(height(transmissionTree.Nodes),1);
    sum(sampleIdx)
    sampledNodes = transmissionTree.Nodes.Name(sampleIdx);
    sampledNodes = sampledNodes(cellfun(@isempty, regexpi(sampledNodes,'ROOT')));
    
tree = samplePhyloTree(transmissionTree,sampledNodes,config);

treeData = exportMetadata(tree,infection, config.files.treeDataFilename);

plot(tree)

%%

R=nan(size(transmissionTree.Nodes{:,1}));

for k = 1:length(transmissionTree.Nodes{:,1})
    R(k)=sum(strcmp(transmissionTree.Edges{:,1}(:,1),transmissionTree.Nodes{k,1}));
    
end
[mean(R(1:50)),mean(R(1:50)==0)]


