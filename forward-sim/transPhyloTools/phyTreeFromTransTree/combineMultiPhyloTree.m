function connectedAlmostPhyloTree = combineMultiPhyloTree(almostPhyloTree)
% take all roots and join by date

% find connected trees within total digraph and their root nodes
[bin,~] = conncomp(almostPhyloTree,'Type','weak');
rootNodes=cell(1,max(bin));
for k=1:max(bin)
    idx = bin == k;
    SG = subgraph(almostPhyloTree, idx);
    node = SG.Nodes.Name{1};
    predNode = predecessors(SG,node);
    while ~isempty(predNode) 
        node = predNode{1};
        predNode = predecessors(SG,node);
    end
    rootNodes{k}=node;
end


connectedAlmostPhyloTree = almostPhyloTree;

% add all roots with weight resembling typical near-neighbor weights
% this is a hacky psuedo-coalsecent model. Should be done right!
% Better hack would be to randomly join tips and internal nodes according
% to typical distances, instead of just joining internal (root) nodes

% find typical pairwise distances between ancestral nodes
nearestDistances=[];
for k=1:max(bin)
    idx = bin == k;
    SG = subgraph(almostPhyloTree, idx);
    
    % find typical nearest neighbor distances
    notSampledNodes=SG.Nodes.Name(outdegree(SG)~=0 & cellfun(@isempty,regexpi(SG.Nodes.Name,'ROOT')));
    
    SG.Edges.Weight(SG.Edges.Weight==0)=1e-4;
    A=adjacency(SG,'weighted');
    A=A+transpose(A);
    uSG = graph(A, SG.Nodes.Name);

    D=distances(uSG,notSampledNodes);
    D(D<2e-4)=1e8;
    
    if min(size(D))>1
        nearestDistances =[nearestDistances,min(D,[],1)];
    end
end
nearestDistances=nearestDistances(nearestDistances~=1e8);

% sample weights to accumulate
if length(rootNodes)> length(nearestDistances)  % very short sim problem
    weight=randsample(nearestDistances/2,length(rootNodes),true);
else
    weight=randsample(nearestDistances/2,length(rootNodes));
end
weight(1)=0;

% iterate through root nodes, adding new node, with sampled weights
node1=rootNodes{1};
for k=2:length(rootNodes)
    node2=rootNodes{k};
    newNode = ['JOIN_',node1,'_',node2];
    
    % handle date differences in rootnewNode=[node,'_anc',num2str(n)];
    node1Date = connectedAlmostPhyloTree.Nodes.timeInfected(ismember(connectedAlmostPhyloTree.Nodes.Name,node1));
    node2Date = connectedAlmostPhyloTree.Nodes.timeInfected(ismember(connectedAlmostPhyloTree.Nodes.Name,node2));
    
    newNodeDate = min([node1Date,node2Date])-weight(k);
    
    newNodeProps = table({newNode}, newNodeDate, 'VariableNames', {'Name' 'timeInfected'});
    connectedAlmostPhyloTree = addnode(connectedAlmostPhyloTree,newNodeProps);
    
    connectedAlmostPhyloTree = addedge(connectedAlmostPhyloTree,newNode,node1,node1Date-newNodeDate);
    connectedAlmostPhyloTree = addedge(connectedAlmostPhyloTree,newNode,node2,node2Date-newNodeDate);
   
    node1 = newNode;    
    
end

end