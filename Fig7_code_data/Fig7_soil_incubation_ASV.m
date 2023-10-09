clear all
close all

rng(1)
load('ASV_dataTable.mat')

%extract relevant data from table
scores = table2array(dataTable(:,contains(dataTable.Properties.VariableNames,'zScore')));
T0_abund = table2array(dataTable(:,contains(dataTable.Properties.VariableNames,'T0')));
soils = {'soil03','soil05','soil06','soil09','soil11','soil12','soil14','soil15','soil16','soil17'};
phylum = dataTable.Phylum;
class = dataTable.Class;
order = dataTable.Order;
genus = dataTable.Genus;
ASV   = dataTable.ASV;
seq   = dataTable.Row;

alpha = 0.01;
enriched = zeros(size(scores));
for i = 1:30
    %minimum values are pseudocounts. Use this information to estimate
    %species richness
    N_asv(i) = sum( (T0_abund(:,i) - min(T0_abund(:,i)))>0 );
    
    %compute zscore threshold using 1) a 1% significance level (two tailed
    %test), 2) a Bonferroni correction, where the correction factor is the
    %number of ASVs, 3) the inverse standard normal cdf
    zthresh(i) = norminv(1-alpha/2/N_asv(i));
    
    %identify which ASVs are positively enriched
    enriched(:,i) = scores(:,i)>zthresh(i);
end

%% Read in results from picrust2 analysis and classify each enriched ASV as a denitrifier, a DNRA strain, or neither
asv_genotype = readtable('denit_DNRA_gene_predicted_count_picrust2.csv','ReadRowNames',1);

for i = 1:height(asv_genotype)
    phy{i} = phylum{find(strcmp(ASV,asv_genotype.ASVID{i}))};
end
asv_genotype.phylum = phy';
asv_genotype.isDenit    = (sum(asv_genotype{:,1:5},2)>0) & (asv_genotype{:,6}~=1);
asv_genotype.isDNRA     = asv_genotype{:,6}==1;
asv_genotype.isNeither  = (~asv_genotype.isDenit) & (~asv_genotype.isDNRA);

%% count the number of ASVs in each phylum and category that are positively enriched
unique_phyla = unique(phylum);

counts_denit = zeros(length(unique_phyla),size(enriched,2));
counts_dnra = zeros(length(unique_phyla),size(enriched,2));
counts_neither = zeros(length(unique_phyla),size(enriched,2));
for i = 1:length(unique_phyla)
    for j = 1:size(enriched,2)
        isDenit = asv_genotype({ASV{enriched(:,j)>0}},:).isDenit;
        isDNRA  = asv_genotype({ASV{enriched(:,j)>0}},:).isDNRA;
        isNeither = asv_genotype({ASV{enriched(:,j)>0}},:).isNeither;
        isPhy = strcmp(asv_genotype({ASV{enriched(:,j)>0}},:).phylum,unique_phyla{i});
        
        counts_denit(i,j) = sum(isDenit & isPhy);
        counts_dnra(i,j) = sum(isDNRA & isPhy);
        counts_neither(i,j) = sum(isNeither & isPhy);
    end
end
counts = counts_denit + counts_dnra + counts_neither;
%%
%sort by the phyla by number of ASVs enriched
sum_counts = sum(counts,2);
[sum_counts,I] = sort(sum_counts,'descend');
counts = counts(I,:);
counts_denit = counts_denit(I,:);
counts_dnra  = counts_dnra(I,:);
counts_neither = counts_neither(I,:);
unique_phyla = unique_phyla(I);

%compute the median across soil replicates
med_nasv           = zeros(1,10);
med_counts         = zeros(size(counts,1),10);
med_counts_denit   = zeros(size(counts,1),10);
med_counts_dnra    = zeros(size(counts,1),10);
med_counts_neither = zeros(size(counts,1),10);
for i = 1:10
    med_nasv(i)             = median( N_asv(((i-1)*3+1):((i-1)*3+3)),2 );
    med_counts(:,i)         = median( counts(:,((i-1)*3+1):((i-1)*3+3)),2 );
    med_counts_denit(:,i)   = median( counts_denit(:,((i-1)*3+1):((i-1)*3+3)),2 );
    med_counts_dnra(:,i)    = median( counts_dnra(:,((i-1)*3+1):((i-1)*3+3)),2 );
    med_counts_neither(:,i) = median( counts_neither(:,((i-1)*3+1):((i-1)*3+3)),2 );
end
sum_med_counts = sum(med_counts,2);

%cut out phyla that are not significantly enriched in any sample
med_counts = med_counts(sum_med_counts>0,:);
med_counts_denit   = med_counts_denit(sum_med_counts>0,:);
med_counts_dnra    = med_counts_dnra(sum_med_counts>0,:);
med_counts_neither = med_counts_neither(sum_med_counts>0,:);
unique_phyla = unique_phyla(sum_med_counts>0);

%% plot
figure()
set(gcf,'units', 'inch','position',[0,0,5,4])
subplot(4,5,[1:4,6:9,11:14])
for i = 1:size(med_counts,1)
    for j = 1:size(med_counts,2)
        if med_counts(i,j) > 0
            scatter(j,i,.00001,0.5*[1,1,1],"filled")
            hold on
            pie_xy([med_counts_denit(i,j),med_counts_dnra(i,j),med_counts_neither(i,j)],j,i,sqrt(med_counts(i,j))*.1135);
        end
    end
end
set(gca,'Ydir','reverse','LineWidth',1,'FontSize',7,'FontName','Arial')
axis image
xlim([0.5,10.5])
ylim([0.5,7.5])
yticks(1:size(med_counts,1))
xticks(1:10)
xticklabels([])
yticklabels(unique_phyla)
ylabel('phylum')
grid on
box on
title('median ASVs enriched per phylum')

%legend for first panel
subplot(4,5,[5,10,15])
dot_weights = 15:-3:0;
for i = 1:6
    if dot_weights(i)>0
        scatter(1,i,dot_weights(i)*25,0.5*[1,1,1],"filled")
    end
    text(1.05,i,num2str(dot_weights(i)))
    hold on
end
set(gca,'Ydir','reverse','YAxisLocation','right','FontSize',7,'FontName','Arial')
ylim([-0.5,6.5])
xlim([0.9 1.1])
xticks([])
yticklabels([dot_weights,0])
title('# ASVs enriched')
%Hide axes lines and ticksq
ax = gca;
ax.Color = 'none';
ax.XAxis.Color = 'none';
ax.YAxis.Color = 'none';
box off;

% bar plots beneath showing fraction of ASVs enriched
subplot(4,5,16:19)
med_frac_asv = sum(med_counts_denit,1)./med_nasv;
bar(1:10,med_frac_asv*100,'FaceColor',[0,0,1],'EdgeColor','none')
for i = 1:10
    text(i,med_frac_asv(i)*100+0.15,[num2str(sum(med_counts_denit(:,i))),'/',num2str(med_nasv(i))],'HorizontalAlignment','center','FontSize',7)
end
set(gca,'LineWidth',1,'FontSize',7,'FontName','Arial')
xlim([0.5,10.5])
ylim([0,1.0])
xticks(1:10)
xticklabels(round(nativePH,1))
ylabel('% ASVs enriched')
xlabel('native soil pH')
ax = gca;
ax.XGrid = 'on';

saveas(gcf,'Fig7BC.svg')