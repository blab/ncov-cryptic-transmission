function phytreeMetadata = exportMetadata(phytreeObj,infectionStruct, filename)
%  output metadata of nodes in phytree

phytreeMetadata.id = get(phytreeObj,'NodeNames');
%cluster
if all(~isempty(regexp(phytreeMetadata.id,'_','split')))
    phytreeMetadata.clusterId=cellfun(@(x) x(1),regexp(phytreeMetadata.id,'_','split'));
end

% real nodes
[ix,loc]=ismember(phytreeMetadata.id,regexprep({infectionStruct.id},':','_'));

names=fieldnames(infectionStruct);
names=names(~ismember(names,{'id','infectedById'}));
for k=1:length(names)
    if isa(infectionStruct(1).(names{k}),'cell')
        phytreeMetadata.(names{k})=nan(size(phytreeMetadata.id));
        phytreeMetadata.(names{k})(ix)=[infectionStruct(loc(ix)).(names{k})]';
    elseif isa(infectionStruct(1).(names{k}),'numeric') || isa(infectionStruct(1).(names{k}),'float') || isa(infectionStruct(1).(names{k}),'logical') 
        phytreeMetadata.(names{k})=nan(size(phytreeMetadata.id));
        phytreeMetadata.(names{k})(ix)=[infectionStruct(loc(ix)).(names{k})]';
    elseif isa(infectionStruct(1).(names{k}),'char')
        phytreeMetadata.(names{k})=cell(size(phytreeMetadata.id));
        phytreeMetadata.(names{k})(ix)={infectionStruct(loc(ix)).(names{k})}';
    end
end

% match inserted nodes with nearest real nodes
for k=find(~ix)'
    if ~isempty(regexpi(phytreeMetadata.id{k},'anc')) && isempty(regexpi(phytreeMetadata.id{k},'JOIN')) && isempty(regexpi(phytreeMetadata.id{k},'ROOT'))  
        % inserted ancestors with 0 distance from true ancestor
        tmpNode = regexprep(phytreeMetadata.id(k),'_anc.+','');
        idx =ismember(regexprep({infectionStruct.id},':','_'),tmpNode);
        for n=1:length(names)
            if isa(infectionStruct(1).(names{n}),'numeric') || isa(infectionStruct(1).(names{n}),'float') || isa(infectionStruct(1).(names{n}),'cell') || isa(infectionStruct(1).(names{n}),'logical')
                phytreeMetadata.(names{n})(k)=infectionStruct(idx).(names{n});
            elseif isa(infectionStruct(1).(names{n}),'char')
                phytreeMetadata.(names{n}){k}=infectionStruct(idx).(names{n});
            end
        end
    elseif ~isempty(regexpi(phytreeMetadata.id{k},'JOIN'))
        % coalescent ancestors from joining pedigrees at start of simulation
    end
    
end

phytreeMetadata = struct2table(phytreeMetadata);

writetable(phytreeMetadata,filename)

end
