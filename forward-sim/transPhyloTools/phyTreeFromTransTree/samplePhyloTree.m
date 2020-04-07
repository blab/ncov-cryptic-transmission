function [timeTree, divTree] = samplePhyloTree(transmissionTree,sampledNodes,config)
% main function of phyTreeFromTransTree for sampling phylogenetic tree from
% transmission tree. 

downsampledTransmissionTree = downsampleTransTree(transmissionTree, sampledNodes);

multiPhyTree = phyloTreeFromTransTree(downsampledTransmissionTree,sampledNodes);

connectedAlmostPhyloTree = combineMultiPhyloTree(multiPhyTree);

treStr = bifurcatingDigraphToNewick(connectedAlmostPhyloTree,config.files.timeTreeFilename);
timeTree=phytreeread(treStr);

end