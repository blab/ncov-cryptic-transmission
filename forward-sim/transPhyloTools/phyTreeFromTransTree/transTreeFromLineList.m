function transmissionTree = transTreeFromLineList(infectionStruct, directed)
% construct transmission tree in nice matlab digraph format (directed graph with edgelist, nodelist)

    if nargin==1
        directed=true;
    end
    
    % convert to unweighted transmission tree graph to get edge list and node list
    infectionGraph=digraph({infectionStruct.infectedById},{infectionStruct.id});

    % add sampling time data to node data
    [ix,loc]=ismember(infectionGraph.Nodes.Name,{infectionStruct.id});
    timeInfected=nan(length(infectionGraph.Nodes.Name),1);
    timeInfected(ix)=[infectionStruct(loc(ix)).timeInfected];
    % give roots time of first infection
    for k=find(~ix)'
        node = infectionGraph.Nodes.Name{k};
        childNode = infectionStruct(ismember({infectionStruct.infectedById},node)).id;
        timeInfected(k) = infectionStruct(ismember({infectionStruct.id},childNode)).timeInfected;
    end
    
    % find distances between nodes as weights
    [~,locParent]=ismember(infectionGraph.Edges.EndNodes(:,1),infectionGraph.Nodes.Name);
    [~,locChild]=ismember(infectionGraph.Edges.EndNodes(:,2),infectionGraph.Nodes.Name);
    weights = timeInfected(locChild) - timeInfected(locParent);
    
    % time-interval weighted directed graph transmission tree
    edgeTable = table(infectionGraph.Edges{:,1},weights, ...
        'VariableNames',{'EndNodes','Weight'});
    nodeTable= table(infectionGraph.Nodes.Name, timeInfected, ...
        'VariableNames',{'Name','timeInfected'});

    if directed
        transmissionTree = digraph(edgeTable,nodeTable);
    else
        transmissionTree = graph(edgeTable,nodeTable);
    end

end