function [status,cmdOut]=prettyRtreePlot(params)
% wrapper to ggtree plot

workDir=cd;
currDir=regexprep(mfilename('fullpath'),'prettyRtreePlot','');
cd(currDir);

cmdStr = 'Rscript --vanilla ggplotPhyloTree.R';    

names=fieldnames(params);
for k=1:length(names)
    if ~isempty(regexp(names{k},'File', 'once'))
        params.(names{k})=[workDir,'/',params.(names{k})];
    end
    cmdStr=[cmdStr,' --',names{k},'="',params.(names{k}),'"'];
end
cmdStr=regexprep(cmdStr,'\\','/');

[status,cmdOut]=system(cmdStr);

cd(workDir);

end
