function divTree = transformTimeTreeToDivergenceTree(timeTree, mutationRate, genomeLength,filename)
% converts time tree to simulated divergence tree with fixed mutation rate

timeDat = get(timeTree);

D=poissrnd(timeDat.Distances*mutationRate*genomeLength);

divTree = phytree(timeDat.Pointers,D,timeDat.NodeNames);

if nargin==4
    phytreewrite(filename,divTree);
    str=fileread(filename);
    str=regexprep(str,'\r\n','');
    fid=fopen(filename,'w+');
    fprintf(fid,str);
    fclose(fid);
end

end