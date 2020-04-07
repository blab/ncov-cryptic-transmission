function almostPhyloTree = phyloTreeFromTransTree(downsampledTransmissionTree,sampledNodes)
% manipulate to get closer to a (bifurcating) phylogenetic tree (with terminal sampled nodes)


almostPhyloTree = downsampledTransmissionTree;

% if sampled node has descendants, convert to outdegree=0 by adding ancestor node and edge with 0 weight
% convert all internal nodes with outdegree>2 into degree 2 by adding fake nodes 
numDescendants=outdegree(almostPhyloTree);
tmp=sum(numDescendants==0);

idx=find(ismember(almostPhyloTree.Nodes.Name,sampledNodes) | numDescendants >2)';
for k=idx
    node = almostPhyloTree.Nodes.Name{k};

    % find descendent nodes
    descNode = successors(almostPhyloTree,node);
    
    n=0; % new node counter
    while ~isempty(descNode)
        n=n+1;
        % find ancestor node at this tree iteration 
        ancNode = predecessors(almostPhyloTree,node);

        % do one node at a time
        descNode=descNode{1};

        % add new ancestor node to descendant with weight equal to original weight
        % and weight to sampled node of 0
        newNode=[node,'_anc',num2str(n)];
        newNodeDate = downsampledTransmissionTree.Nodes.timeInfected(ismember(downsampledTransmissionTree.Nodes.Name,node));
        newNodeProps = table({newNode}, newNodeDate, 'VariableNames', {'Name' 'timeInfected'});
        almostPhyloTree = addnode(almostPhyloTree,newNodeProps);

        wIdx=ismember(almostPhyloTree.Edges{:,1}(:,1),node) & ismember(almostPhyloTree.Edges{:,1}(:,2),descNode);
        w=almostPhyloTree.Edges.Weight(wIdx);
        almostPhyloTree = addedge(almostPhyloTree,newNode,descNode,w);
        almostPhyloTree = addedge(almostPhyloTree,newNode,node,0);

        % remove edges between original sampled node and descendant
        almostPhyloTree = rmedge(almostPhyloTree,node,descNode);

        if ~isempty(ancNode)  % if node isn't root
            % add links from grandparent node to new node with original weight
            wIdx=ismember(almostPhyloTree.Edges{:,1}(:,1),ancNode) & ismember(almostPhyloTree.Edges{:,1}(:,2),node);
            w=almostPhyloTree.Edges.Weight(wIdx);
            almostPhyloTree = addedge(almostPhyloTree,ancNode,newNode,w); 

            % remove edges between original sampled node and grandparent
            almostPhyloTree = rmedge(almostPhyloTree,ancNode,node);
        end
        
        % repeat for other descendents 
        descNode = successors(almostPhyloTree,node);
    end
    display(['make bifurcating: percent complete: ',num2str(100*k/max(idx))])
end

% remove residual unsampled tips with outdegree 0 (these are descendants of multifurcating nodes handled earlier)
    numDescendants=outdegree(almostPhyloTree);

    % find terminal nodes that aren't sampled
    remNodes = setdiff(almostPhyloTree.Nodes.Name(numDescendants==0),sampledNodes);
    
    % find edges that end in unsampled terminal nodes
    remEdges=find(ismember(almostPhyloTree.Edges{:,1}(:,2), remNodes));

    almostPhyloTree = rmedge(almostPhyloTree,remEdges);
    almostPhyloTree = rmnode(almostPhyloTree,remNodes);

    
% a phylogenetic tree only has nodes of outdegree==0 (tips) and outdegree == 2.
% aggregate over internal nodes with outdegree == 1 and give total weight
% to outdegree 2 ancestors
numDescendants=outdegree(almostPhyloTree);
idx=find(ismember(almostPhyloTree.Nodes.Name,sampledNodes) | numDescendants ==2)';
for k=idx
    node = almostPhyloTree.Nodes.Name{k};
    ancNode = predecessors(almostPhyloTree,node);
    
    while ~isempty(ancNode) 
        descAncNode = successors(almostPhyloTree,ancNode);

        if length(descAncNode)==1
            predAncNode = predecessors(almostPhyloTree,ancNode);

            if ~isempty(predAncNode)
                edges=table2array(almostPhyloTree.Edges(:,1));
                wIdx=ismember(edges(:,1),predAncNode) & ismember(edges(:,2),ancNode);
                w=almostPhyloTree.Edges.Weight(wIdx);

                wIdx=ismember(edges(:,1),ancNode) & ismember(edges(:,2),node);
                w = w + almostPhyloTree.Edges.Weight(wIdx);

                almostPhyloTree = addedge(almostPhyloTree,predAncNode,node,w);

                almostPhyloTree = rmedge(almostPhyloTree,predAncNode,ancNode);
                almostPhyloTree = rmedge(almostPhyloTree,ancNode,node);
            end
            ancNode=predAncNode;
        else
            ancNode={};
        end
    end
    display(['sum branches: percent complete: ',num2str(100*k/max(idx))])    
end
% get rid of nodes that aren't in edges
almostPhyloTree = rmnode(almostPhyloTree,setdiff(almostPhyloTree.Nodes.Name,unique(almostPhyloTree.Edges{:,1})));

% get rid of remaing nodes with outdegree == 1 (roots that have no descendents)
roots = almostPhyloTree.Nodes.Name(outdegree(almostPhyloTree)==1);
rmEdgeIdx = find(ismember(almostPhyloTree.Edges{:,1}(:,1),roots)); 

almostPhyloTree = rmedge(almostPhyloTree, rmEdgeIdx);
almostPhyloTree = rmnode(almostPhyloTree, roots);



end