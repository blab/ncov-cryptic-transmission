
clear; clc;

% include
addpath(genpath('../transPhyloTools'));

%% today's date time span

datenum('15-Mar-2020')/365.2431 - 2020.038161 


%% run model
config=yaml.ReadYaml('configTreeClusterSim_snohomish_1.yaml');

T=config.simulation.startTime:config.simulation.timeStep:(config.simulation.startTime+ config.simulation.totalTime);
Tstr=datestr(T*365.2431,'dd-mmm');
x_ticks = find(ismember(Tstr, {'15-Jan','01-Feb','15-Feb','01-Mar','15-Mar'}));


numReps=500;


finalIncidence=nan(numReps,1);
finalPrevalence = nan(numReps,1);
doublingTime = nan(numReps,1);
compScore = nan(numReps,1);

incidence=nan(numReps,length(T));
prevalence=nan(numReps,length(T));

nClust=nan(numReps,1);
nBigClust=nan(numReps,1);

count=0;
for k=1:numReps
    
    endPrev=0;
    while endPrev<7 
        totalInfections=treeClusterSim(config);
        
        infection=xls2struct(config.files.infectionLinelistFilename,'structArray');
        clust=xls2struct(config.files.clusterLinelistFilename,'structArray');

        [incidence(k,:),x]=hist([infection.timeInfected],T);
        incidence(k,:) = cumsum(incidence(k,:));
        finalIncidence(k)=incidence(k,end);
        
        recovered = cumsum(hist([infection.timeInfected] + [infection.latentDuration] + [infection.infectiousDuration],[T,T(end)+config.simulation.timeStep]));
        prevalence(k,:) = incidence(k,:)-recovered(1:end-1);

        mod = fitlm(x(end-20:end),log2(incidence(k,end-20:end)));
        doublingTime(k)=365.2431/mod.Coefficients.Estimate(2);

        finalPrevalence(k)=prevalence(k,end);
        endPrev =  prevalence(k,end-31);
      
        if endPrev==0
            count=count+1;
        end
        
        
        disp('done!')


        nClust(k)=length(clust);
        nBigClust(k)=sum([clust.totalInfections]>=10);

        
       
    end
    

end


%% paper figure

clc
cut_dates = find(ismember(Tstr, {'01-Mar','05-Mar','10-Mar','15-Mar'}));

figure(20); clf;

str = '#4141DD';
str = '#666666';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

medPrev=quantile(prevalence(:,end),0.50);
prevScore = (prevalence(:,end)-medPrev).^2;
minCompIdx=find(prevScore==min(prevScore),1);

subplot(1,2,1); hold on;
plot(T,prevalence,'color',color);
ylabel('prevalence')
% set(gca,'yscale','log');
plot(T,prevalence(minCompIdx,:),'r','linewidth',1.5);
set(gca,'xtick',T(x_ticks), 'xticklabel',Tstr(x_ticks,:))
axis tight
box off
ax=get(gca,'ylim');

today = find(ismember(Tstr, {'05-Mar'}));

plot(T(today)*[1,1],ax,'k--')



% subplot(1,2,2); hold on;
% hist(finalPrevalence,20);
% y=hist(finalPrevalence,20,'facecolor',color,'color','none');
% xlabel('prevalence')
% title(['mean = ',num2str(mean(finalPrevalence))], 'fontweight','normal')

x=0:200:8000;
for k=4:-1:1
    subplot(8,2,4*k); hold on;
    [y,x]=hist(prevalence(:,cut_dates(k)),x);
    bar(x,y,1,'facecolor',color,'edgecolor','none' );
    round([nanmedian(prevalence(:,cut_dates(k))), quantile(prevalence(:,cut_dates(k)),.05), quantile(prevalence(:,cut_dates(k)),0.95)])
    round([nanmedian(incidence(:,cut_dates(k))), quantile(incidence(:,cut_dates(k)),.05), quantile(incidence(:,cut_dates(k)),0.95)])
    title(Tstr(cut_dates(k),:), 'fontweight','normal')
    box off
    axis tight
    if k==4
        ax=get(gca,'xlim');
        ax(1)=0;
        xlabel('prevalence')

    end
    set(gca,'xlim',ax,'ytick',[]);
    
end



figure(21); clf;
subplot(1,2,1);

x=min(nClust):max(nClust);
y=hist(nClust,x);
bar(x,y,1,'facecolor',color,'edgecolor','none' );
box off
axis tight
    
subplot(1,2,2);

x=min(nBigClust):max(nBigClust);
y=hist(nBigClust,x);
bar(x,y,1,'facecolor',color,'edgecolor','none' );
box off
axis tight
([nanmean(nClust-1), quantile(nClust-1,.05), quantile(nClust-1,0.95)])
([nanmean(nBigClust-1), quantile(nBigClust-1,.05), quantile(nBigClust-1,0.95)])


save('paper_run.mat','prevalence','incidence')
