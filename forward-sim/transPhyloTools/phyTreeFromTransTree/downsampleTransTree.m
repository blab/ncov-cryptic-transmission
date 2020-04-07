function transmissionTreeObj = downsampleTransTree(transmissionTreeObj, sampledNodes)
% downsample tree to only include ancestry of sampled nodes

% identify all infections in sample ancestry by traversing tree backward from samples
inSampleIdx = false(height(transmissionTreeObj.Nodes),1);
for k=1:length(sampledNodes)
    ancestralNode = sampledNodes{k};
    while ~isempty(ancestralNode) && ~ismember(ancestralNode,transmissionTreeObj.Nodes.Name(inSampleIdx))
        inSampleIdx(ismember(transmissionTreeObj.Nodes.Name,ancestralNode))=true;
        ancestralNode = predecessors(transmissionTreeObj,ancestralNode);
    end
    display(['downsampling: percent complete: ',num2str(100*k/length(sampledNodes))])
end

rmEdgeIdx = ~ismember(transmissionTreeObj.Edges{:,1}(:,1),transmissionTreeObj.Nodes.Name(inSampleIdx)) ...
    | ~ismember(transmissionTreeObj.Edges{:,1}(:,2),transmissionTreeObj.Nodes.Name(inSampleIdx));

rmEdgeIdx = find(rmEdgeIdx);

transmissionTreeObj = rmedge(transmissionTreeObj,rmEdgeIdx);

rmNodes = setdiff(transmissionTreeObj.Nodes.Name,transmissionTreeObj.Nodes.Name(inSampleIdx));
transmissionTreeObj = rmnode(transmissionTreeObj,rmNodes);

end


