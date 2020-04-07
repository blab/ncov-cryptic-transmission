function [infection, cluster]=plotTreeClusterSim(config)

if isstruct(config)
    
elseif ischar(config)
    config=yaml.ReadYaml(config);
else 
    error('config must be struct or filename');
end

cluster=xls2struct(config.files.clusterLinelistFilename,'structArray');
infection=xls2struct(config.files.infectionLinelistFilename,'structArray');

%% clusters



[~,sortIdx]=sort([cluster.timeCreated]);
cluster=cluster(sortIdx);

[uniqueClusters, ~,clusterOrderIdx]=unique({cluster.id});
for k=1:length(cluster)
    cluster(k).smallId=clusterOrderIdx(k);
    cluster(k).smallSeededById=find(ismember(uniqueClusters,cluster(k).seededById));
end



figure(1); clf

subplot(2,2,1); hold all;
for k=1:length(cluster)
    plot([cluster(k).timeCreated,cluster(k).timeKilled],k*[1,1])
    ylabel('cluster');
    title('Cluster duration','fontweight','normal');
end

subplot(2,2,2); hold all;
[y,x]=hist(log10([cluster.totalInfections]),linspace(0,log10(max([cluster.totalInfections])),10));
bar(x,y); 
xlabel('log10 total infections'); ylabel('number of clusters')
title('SIR cluster size','fontweight','normal')

subplot(2,2,3);
transIdx=cellfun(@isempty,regexp({cluster.seededById},'ROOT'));
A=sparse([cluster(transIdx).smallSeededById],[cluster(transIdx).smallId] ...
    ,ones(sum(transIdx),1),length(cluster),length(cluster));

clusterGraph=digraph({cluster(transIdx).seededById},{cluster(transIdx).id});
plot(clusterGraph)
title('Importations: connected clusters','fontweight','normal')




%% individual

[y,x]=hist([infection.timeInfected],config.simulation.startTime+(0:config.simulation.timeStep:config.simulation.totalTime));

prevalence=0*x;
for k=1:length(x)
    prevalence(k)=sum([infection.timeInfected]<=x(k)) - sum(([infection.timeInfected] + [infection.latentDuration] + [infection.infectiousDuration]) <=x(k));
end


subplot(2,2,4); cla; hold on;
plot(x,cumsum(y),'b');
plot(x,prevalence,'r')
ylabel('count per day')
title('True incidence','fontweight','normal')
mod = fitlm(x,log2(cumsum(y)));

disp(['doubling time = ',num2str(365.2431/mod.Coefficients.Estimate(2)),' days'])

drawnow;


%% individual


figure(2) 
prevalence=0*x;
for k=1:length(x)
    prevalence(k)=sum([infection.timeInfected]<=x(k)) - sum(([infection.timeInfected] + [infection.latentDuration] + [infection.infectiousDuration]) <=x(k));
end


subplot(2,2,1); cla; hold on;

plot(x,prevalence,'r')
ylabel('count per day')
title('Prevalance','fontweight','normal')
mod = fitlm(x,log2(cumsum(y)));

subplot(2,2,2); cla; hold on;
% 
% plot(x,prevalence*400,'r')
% ylabel('count per day')
% title('Population exposed','fontweight','normal')
% mod = fitlm(x,log2(cumsum(y)));


disp(['doubling time = ',num2str(365.2431/mod.Coefficients.Estimate(2)),' days'])

drawnow;


%% save
figure(1);
if isfield(config.files,'incidenceSummaryFilename')
    print(config.files.incidenceSummaryFilename,'-dpng','-r300')
end

end