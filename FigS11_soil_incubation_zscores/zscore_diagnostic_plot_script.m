% Only for dpH=0, for the 10 soils where we have sequencing data.  Can you 
% give us a list of the taxonomy of responding ASVs where response is 
% defined for nitrate+ vs nitrate-. Namely, we want those ASVs that 
% we can say definitely are responding to the presence of nitrate in 
% unperturbed soils. 

clear
close all

% Full dataset filtered to the relevant dpH=0
load('dataSubset_over1000_spikeRemoved_dpH0.mat','d','m','spikeIn');

cts = d{:,1:size(m,1)}+0.5;
cts = cts./spikeIn;
%%
[soilNum, ~, soilIdx] = unique(m.SoilNum);

timepoint = cellfun(@(s)double(s(end)-'0'), m.Time_point);

% For each soil: 
% * 3 replicates of T0
% * 3 replicates with no nitrate (T9)
% * 3 replicates with nitrate (T9)
% * 3 zScores from the 3 nitrate / no-nitrate comparisons (T9)
dataTable = [repmat(d(:,1), [1, 12*max(soilIdx)]), d(:, size(m,1)+1:end)];

% Soils 1..18
native_soil_pH = [4.703, 5.094, 4.987, 5.277, 5.324, 5.405, 5.514, 5.552,...
        5.822, 5.975, 6.186, 6.255, 6.435, 6.545, 6.789, 6.860, 7.052, 7.681];
nativePH = native_soil_pH(soilNum);

SoilID = {'CE239','CE201','CE73','CE277','CE234','CE229','Neutral2','Neutral5','Neutral6','Neutral3'};

cFrac = NaN(1,max(soilIdx));
c0 = NaN(1,max(soilIdx));
for ss=1:max(soilIdx)
    thisSoil = soilIdx==ss;
    thisSoilNum = soilNum(ss);

    % Find samples from this soil with and wihout nitrate
    withNitrate = thisSoil & m.Nitrate_input>0 & timepoint==9;
    noNitrate = thisSoil & m.Nitrate_input==0 & timepoint==9;
    T0 = thisSoil & timepoint==0;
    ctsNitrate = cts(:,withNitrate);
    ctsCtrl = cts(:,noNitrate);
    ctsT0 = cts(:,T0);
    
    % Learn error model by comparing control to itself (between replciates)
    [~, cFracNull(ss), c0Null(ss)] = scoreOutliers(ctsCtrl,ctsCtrl(:,[2 3 1]));
    [~, cFracNO3(ss), c0NO3(ss)] = scoreOutliers(ctsNitrate,ctsNitrate(:,[2 3 1]));
    [~, cFrac(ss), c0(ss)] = scoreOutliers([ctsCtrl,ctsNitrate], [ctsCtrl(:,[2 3 1]), ctsNitrate(:,[2 3 1])]);
    
    % Use these parameters when comparing to nitrate
    z = scoreOutliers(ctsCtrl,ctsNitrate, cFrac(ss), c0(ss));

    %
    varnamesNitrate = compose(sprintf('Soil%d_T9abd_withNitrate_rep%%d',thisSoilNum), 1:3);
    varnamesCtrl = compose(sprintf('Soil%d_T9abd_noNitrate_rep%%d',thisSoilNum), 1:3);
    varnamesT0 = compose(sprintf('Soil%d_T0abd_rep%%d',thisSoilNum), 1:3);
    varnamesZ = compose(sprintf('Soil%d_T9comparison_zScore_rep%%d',thisSoilNum), 1:3);
    idx = 12*(ss-1)+(1:3);

    dataTable.Properties.VariableNames(idx) = varnamesT0;
    dataTable{:,idx} = ctsT0;
    idx = idx+3;
    dataTable.Properties.VariableNames(idx) = varnamesCtrl;
    dataTable{:,idx} = ctsCtrl;
    idx = idx+3;
    dataTable.Properties.VariableNames(idx) = varnamesNitrate;
    dataTable{:,idx} = ctsNitrate;
    idx = idx+3;
    dataTable.Properties.VariableNames(idx) = varnamesZ;
    dataTable{:,idx} = z;

    %
    if ss < 10
        ax = subplot(4,3,ss);
    else
        ax = subplot(4,3,11);
    end
    rg = [5e-5,10];
    col = min(abs(z),5);
    plot(rg,rg,'k-');
    hold on
    colormap jet
    scatter(ctsCtrl(:),ctsNitrate(:),10,col(:),'filled')
    ax.XScale = 'log';
    ax.YScale = 'log';
    axis([rg rg])
    axis square
    title({sprintf([SoilID{ss},' (pH = %.1f)'],nativePH(ss)),sprintf('Error model = %.0f%% + %.1f',100*cFrac(ss),c0(ss))});
    xlabel('abundance (control)');
    ylabel('abundance (nitrate^+)');
    h=colorbar();
    h.Label.String = "z-score";
    
end
fprintf('** DIAGNOSTICS **\n');
fprintf('          (nitrate)      (no nitrate)      (combined)\n');
fprintf('cFrac:   %.2f +- %.2f    %.2f +- %.2f    %.2f +- %.2f\n', mean(cFracNO3), std(cFracNO3), mean(cFracNull), std(cFracNull), mean(cFrac), std(cFrac));
fprintf('c0:      %.2f +- %.2f    %.2f +- %.2f    %.2f +- %.2f\n', mean(c0NO3), std(c0NO3), mean(c0Null), std(c0Null), mean(c0), std(c0));
set(gcf,'Position',[100 100 800 960])
%%
print(gcf,'FigSX_zscore.png','-dpng','-r150');
save('ASV_dataTable.mat', 'dataTable', 'nativePH')
writetable(dataTable,'ASV_dataTable.csv','WriteRowNames',true)