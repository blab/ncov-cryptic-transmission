% builds directed graph equivalent to phylogenetic tree from a linelist
% with who-infected-whom data.
% then move to R to visualize. 

clc; clear; 

addpath(genpath('../utils'))

%% load example line list
% must have fields "id","infectedById","timeInfected"

infection=xls2struct(['example/infectionLineList.csv'],'structArray');

%% build transmission tree from line list
transmissionTree = transTreeFromLineList(infection)

%% select nodes to observe and downsample tree to only include ancestry of observed

% sample nodes
rng(1000)
sampleIdx=(1/200*ones(height(transmissionTree.Nodes),1))>rand(height(transmissionTree.Nodes),1);
sum(sampleIdx)
sampledNodes = transmissionTree.Nodes.Name(sampleIdx);

downsampledTransmissionTree = downsampleTransTree(transmissionTree, sampledNodes)

%% convert to (possibly set of disconnected) bifurcating trees with only sampled nodes at tips
multiPhyTree = phyloTreeFromTransTree(downsampledTransmissionTree,sampledNodes)

%% join up components into one big tree

connectedAlmostPhyloTree = combineMultiPhyloTree(multiPhyTree)


% this is a simulated complete phylogenetic tree! But it a digraph object
% and not a phytree object...

%% write to newick and get as phytree

str = bifurcatingDigraphToNewick(connectedAlmostPhyloTree,'example/connectedTree.tre');

tre=phytreeread('example/connectedTree.tre');
plot(tre)

%% export metadata
treMetadata = exportMetadata(tre,infection, 'example/connectedTree.csv');

%% export pretty tree from R
params.treeFile = 'example/connectedTree.tre';
params.dataFile = 'example/connectedTree.csv';
params.legendPosition = 'right';

params.colorField = 'clusterId';
prettyRtreePlot(params)

params.colorField = 'timeInfected';
prettyRtreePlot(params)


%% create fake divergence tree

divTree = transformTimeTreeToDivergenceTree(tre, 0.01,906,'example/connectedTree.divergence.tre' );

params.treeFile = 'example/connectedTree.divergence.tre';
params.dataFile = 'example/connectedTree.csv';
params.legendPosition = 'right';
params.scale = 'divergence';

params.colorField = 'timeInfected';
prettyRtreePlot(params)
















