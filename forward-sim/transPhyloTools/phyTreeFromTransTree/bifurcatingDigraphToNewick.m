function [treeStr, subtreeStr] = bifurcatingDigraphToNewick(treeGraph,filename)

% check if valid bifurcating tree
    numDescendents = outdegree(treeGraph);
    if any(numDescendents~=0 & numDescendents~=2) && numDescendents(end)~=1
        error('invalid outdegree','must be bifurcating tree (possibly with singleton root).')
    else
        numLeaves = sum(numDescendents==0);
        numNodes  = sum(numDescendents==2);
        if (numLeaves - numNodes) ~= 1
            error('invalid tree','must have numLeaves = numNodes+1')
        end
    end
    
% format names with branch length
    names = treeGraph.Nodes.Name;

    [ix,loc]=ismember(names,treeGraph.Edges{:,1}(:,2));
    
    % add replace priveleged characters in the Newick format
    names = regexprep(names,':','_');
    names = regexprep(names,',','-');
    names = regexprep(names,';','+');
    
    branchLengths(ix) =  treeGraph.Edges{:,2}(loc(ix));
    branchLengths(~ix) = 0;
    for k=1:length(names)
        names{k}=[names{k},':',num2str(branchLengths(k))];
    end
    
% construct node subtree strings
    subtreeStr  = cell(numLeaves+numNodes,1);
    
% construct leaf strings
    subtreeStr(numDescendents==0)=names(numDescendents==0);

% find depth of nodes for breadth-first traversal
    depth=nan(treeGraph.numnodes,1);
    for k=1:treeGraph.numnodes
        count=0;
        node = treeGraph.Nodes.Name{k};
        ancNode = predecessors(treeGraph,node);
        while ~isempty(ancNode)
            ancNode = predecessors(treeGraph,ancNode);
            count=count+1;
        end
        depth(k)=count;
    end

% order internal nodes, starting from deepest
    nodes = treeGraph.Nodes.Name(numDescendents>0);
    [~,nodeOrder] = sort(depth(numDescendents>0),'descend');
    nodes = nodes(nodeOrder);

% traverse tree and build Newick subtree strings
    for k = 1:length(nodes)
        children = successors(treeGraph,nodes{k});
        for n=1:length(children)
            childNodeIdx=ismember(treeGraph.Nodes.Name,children(n));
            if numDescendents(childNodeIdx)==2
                grandChildren = successors(treeGraph,children{n});
                grandChildrenIdx = find(ismember(treeGraph.Nodes.Name,grandChildren));
                grandChildStr = ['(',subtreeStr{grandChildrenIdx(1)},',',subtreeStr{grandChildrenIdx(2)},')'];
                
                subtreeStr{childNodeIdx} = [grandChildStr,names{childNodeIdx}];
             end
        end
    end
    % root
    childrenIdx = find(ismember(treeGraph.Nodes.Name,children));
    for k=1:length(names)
        if isempty(regexp(subtreeStr{childrenIdx(1)},names{k}, 'once')) && isempty(regexp(subtreeStr{childrenIdx(2)},names{k}, 'once'))
            rootNode=names{k};
        end
    end
    rootSubtreeIdx = find(cellfun(@isempty,subtreeStr));
    
    subtreeStr{rootSubtreeIdx} = ['(',subtreeStr{childrenIdx(1)},',',subtreeStr{childrenIdx(2)},')',rootNode,';'];
    
% output    
    treeStr = subtreeStr{rootSubtreeIdx};
    for k=1:length(subtreeStr)-1
        subtreeStr{k}=[subtreeStr{k},';'];
    end


% write to file
if nargin>1
    fid = fopen(filename,'wt') ;
    fprintf(fid,'%s',treeStr) ;
    fclose(fid);
end

end