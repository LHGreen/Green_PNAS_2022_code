% load data

load('PopVec.mat')

dx = 40; % bin size 40 mm

xBins = -200:dx:2000; %0:dx:1800;
xVals = (xBins(2:end) + xBins(1:end-1))/2;

colors = cbrewer('qual', 'Paired', 6);
%% Smooth the data a bit
b = 5;
for ii = 1:10
    for jj = 1:size(CA1spikes{ii}, 1)
        CA1spikes{ii}(jj,:) = movmean(CA1spikes{ii}(jj,:), b, 'omitnan');
    end
end

for ii = 1:10
    for jj = 1:size(CA3spikes{ii}, 1)
        CA3spikes{ii}(jj,:) = movmean(CA3spikes{ii}(jj,:), b, 'omitnan');
    end
end

for ii = 1:10
    for jj = 1:size(LSspikes{ii}, 1)
        LSspikes{ii}(jj,:) = movmean(LSspikes{ii}(jj,:), b, 'omitnan');
    end
end

%% Plot each condition. These are the population vector correlations

figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ii = 1:10
    subplot(3, 4, ii)
imagesc(sort_cells(zscore_nan(CA1spikes{ii}')'))
box off
title(num2str(ii))
end

figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ii = 1:10
    subplot(3, 4, ii)
imagesc(sort_cells(zscore_nan(CA3spikes{ii}')'))
box off
title(num2str(ii))
end

figure('units', 'normalized', 'outerposition', [0 0 1 1])
for ii = 1:10
    subplot(3, 4, ii)
imagesc(sort_cells(zscore_nan(LSspikes{ii}')'))
box off
title(num2str(ii))
end

%% For each region, get indices of cells that have activity greater than background in at least 1 env.
binIDs = 7:length(xBins)-7;%3:size(CA1spikes{1}, 2)-3; % first two bins on either side of the track don't always have data 
noise = 0.1;  % Hz, minimum mean firing rate to include cell
limit = 40; % maximum mean firing rate to include cell, 
% (This includes all cells, E and I, it didn't make much of a
% difference to exclude interneurons)
CA1cells = [];
for ii = 1:10
   CA1cells = [CA1cells find(nanmean(CA1spikes{ii}(:, binIDs)') > noise & nanmean(CA1spikes{ii}(:, binIDs)') < limit)];
    
end
CA1cells = unique(CA1cells);

CA3cells = [];
for ii = 1:10
   CA3cells = [CA3cells find(nanmean(CA3spikes{ii}(:, binIDs)') > noise & nanmean(CA3spikes{ii}(:, binIDs)') < limit)];
    
end
CA3cells = unique(CA3cells);

LScells = [];
for ii = 1:10
   LScells = [LScells find(nanmean(LSspikes{ii}(:, binIDs)') > noise)];
    
end
LScells = unique(LScells);

%% Plot population vectors again with active cells only
f = figure;
f.Position = [10 10 800 800];
for ii = 1:10
    subplot(3, 4, ii)
imagesc(sort_cells(zscore_nan(CA1spikes{ii}(CA1cells, :)')'))
box off
% axis square
title(num2str(ii))
clim([-2 2])
end

f = figure;
f.Position = [10 10 800 800];
for ii = 1:10
    subplot(3, 4, ii)
imagesc(sort_cells(zscore_nan(CA3spikes{ii}(CA3cells, :)')'))
box off
title(num2str(ii))
clim([-2 2])
end

f = figure;
f.Position = [10 10 800 800];
for ii = 1:10
    subplot(3, 4, ii)
    imagesc(sort_cells(zscore_nan(LSspikes{ii}(LScells,:)')'))
    box off
    title(num2str(ii))
    clim([-2 2])
end



%% compare firing rates between conditions

medianFr = nan(10, 1);
meanFr = nan(10, 1);
stdFr = nan(10, 1);
data = [];
for ii = 1:10
    medianFr(ii) = nanmedian(CA1spikes{ii}(:));
    meanFr(ii) = nanmean(CA1spikes{ii}(:));
    stdFr(ii) = nanstd(CA1spikes{ii}(:));
    temp = nanmean(CA1spikes{ii}(CA1cells, binIDs)');
    data(:, ii) = temp(:);
end
% figure
% scatter(1:10, medianFr, 20, 'b', 'filled')
% box off
% xlabel('Condition')
% ylabel('MedianFr')
% 
% figure
% scatter(1:10, meanFr, 20, 'b', 'filled')
% hold on
% scatter(1:10, meanFr + stdFr, 20, 'k', 'filled')


[p, ~, stats] = kruskalwallis(data);
multcompare(stats);

%%
data(isnan(data)) = 0;
figure
hold on
for ii = 1:10
    plot(sort(data(:, ii)))
end
box off
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')

%% Make correlation matrices
% save the corr mats
CA1corr = nan(size(CA1spikes{1}(:, binIDs), 2), size(CA1spikes{1}(:, binIDs), 2), 100);
CA3corr = nan(size(CA3spikes{1}(:, binIDs), 2), size(CA3spikes{1}(:, binIDs), 2), 100);
LScorr = nan(size(LSspikes{1}(:, binIDs), 2), size(LSspikes{1}(:, binIDs), 2), 100);

fig = figure;
fig.Position = [10 10 1200 1200];
for ii = 1:10
    for jj = 1:10
        subplot(10, 10, sub2ind([10 10], ii, jj))
        CA1corr(:,:,sub2ind([10 10], ii, jj)) = corr(zscore(CA1spikes{ii}(CA1cells,binIDs)')', zscore(CA1spikes{jj}(CA1cells,binIDs)')', 'rows', 'complete');
        imagesc( CA1corr(:,:,sub2ind([10 10], ii, jj)))
        box off
        axis square
        xticks([16 31]); xticklabels([60 120])
        yticks([16 31]); yticklabels([60 120])
%         set(gca,'visible','off')
%         set(gca,'xtick',[])
%         colorbar
        clim([-0.5 1])
%         title([num2str(ii) ' ' num2str(jj)])
    end
end


% CA3 
fig = figure;
fig.Position = [10 10 1200 1200];
for ii = 1:10
    for jj = 1:10
        subplot(10, 10, sub2ind([10 10], ii, jj))
        CA3corr(:,:,sub2ind([10 10], ii, jj)) = corr(zscore(CA3spikes{ii}(CA3cells,binIDs)')', zscore(CA3spikes{jj}(CA3cells,binIDs)')', 'rows', 'complete');
        imagesc( CA3corr(:,:,sub2ind([10 10], ii, jj)))
        box off
        axis square
        xticks([16 31]); xticklabels([60 120])
        yticks([16 31]); yticklabels([60 120])
        %         set(gca,'visible','off')
        %         set(gca,'xtick',[])
        %         colorbar
        clim([-0.5 1])
        %         title([num2str(ii) ' ' num2str(jj)])
    end
end

% LS
fig = figure;
fig.Position = [10 10 1200 1200];
for ii = 1:10
    for jj = 1:10
        subplot(10, 10, sub2ind([10 10], ii, jj))
        LScorr(:,:,sub2ind([10 10], ii, jj)) = corr(zscore(LSspikes{ii}(LScells,binIDs)')', zscore(LSspikes{jj}(LScells,binIDs)')', 'rows', 'complete');
        imagesc(LScorr(:,:,sub2ind([10 10], ii, jj)))
        box off
        axis square
         xticks([16 31]); xticklabels([60 120])
        yticks([16 31]); yticklabels([60 120])
%         set(gca,'visible','off')
%         set(gca,'xtick',[])
%         colorbar
        clim([-0.5 1])
%         title([num2str(ii) ' ' num2str(jj)])
    end
end



%% Compare jump trials with run trials before and after
histBins = -1:0.05:1;


temp = [CA1corr(:,:,sub2ind([10 10], 1, 7)); CA1corr(:,:,sub2ind([10 10], 3, 7)); CA1corr(:,:,sub2ind([10 10], 5, 7)); ...
    CA1corr(:,:,sub2ind([10 10], 2, 8)); CA1corr(:,:,sub2ind([10 10], 4, 8)); CA1corr(:,:,sub2ind([10 10], 6, 8))];

temp2 = [CA1corr(:,:,sub2ind([10 10], 1, 9)); CA1corr(:,:,sub2ind([10 10], 3, 9)); CA1corr(:,:,sub2ind([10 10], 5, 9)); ...
    CA1corr(:,:,sub2ind([10 10], 2, 10)); CA1corr(:,:,sub2ind([10 10], 4, 10)); CA1corr(:,:,sub2ind([10 10], 6, 10))];
data = [temp(:); temp2(:)];
numel = length(temp(:));
group = [ones(numel, 1); 2*ones(numel, 1)];

[p, tbl, stats] = kruskalwallis(data, group);

figure
hold on
histogram(temp, histBins, 'Normalization', 'probability', 'FaceColor', colors(1,:))
histogram(temp2, histBins, 'Normalization', 'probability', 'FaceColor', colors(2,:))
box off
xlabel('Correlation Coefficient')
ylabel('Frequency')
title(['p = ' num2str(p)])
legend('Before', 'After')

%%
[N1, ~] = histcounts(temp(:), histBins, 'Normalization', 'Probability');

[N2, ~] = histcounts(temp2(:), histBins, 'Normalization', 'Probability');

fit1 = fitdist(temp(:), 'Kernel');
fit2 = fitdist(temp2(:), 'Kernel');
% xVals = -0.5:0.01:0.8;
xVals = (histBins(2:end) + histBins(1:end-1))/2;
y1 = pdf(fit1, xVals);
y2 = pdf(fit2, xVals);


fig = figure;
subplot(4, 3, 1)
fig.Position = [10 10 800 800];
hold on
histogram(temp, histBins, 'Normalization', 'Probability', 'FaceColor', colors(2,:))
histogram(temp2, histBins, 'Normalization', 'Probability','FaceColor', colors(4,:))
plot(xVals, y1/sum(y1), 'color', colors(2,:), 'lineWidth', 2)
plot(xVals, y2/sum(y2), 'color', colors(4,:), 'lineWidth', 2)
box off
xlabel('Correlation Coefficient')
ylabel('Frequency')
title(['p = ' num2str(p)])
legend('Before', 'After')

% add a QQ plot
figure(fig)
subplot(4, 3, [4 7])
h = qqplot(temp(:));
hold on
h2 = qqplot(temp2(:));

set(h(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(1,:));
set(h(2), 'lineWidth', 2, 'color', colors(2,:))
set(h(3), 'lineWidth', 2, 'color', colors(2,:))
set(h2(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(3,:));
set(h2(2), 'lineWidth', 2, 'color', colors(4,:))
set(h2(3), 'lineWidth', 2, 'color', colors(4,:))


%% Compare correlations in same running direction vs. opposite
% CA1
temp = [CA1corr(:,:,sub2ind([10 10], 1, 3)); CA1corr(:,:,sub2ind([10 10], 1, 5)); CA1corr(:,:,sub2ind([10 10], 3, 5)); ...
    CA1corr(:,:,sub2ind([10 10], 2, 4)); CA1corr(:,:,sub2ind([10 10], 2, 6)); CA1corr(:,:,sub2ind([10 10], 4, 6))];

temp2 = [CA1corr(:,:,sub2ind([10 10], 1, 2)); CA1corr(:,:,sub2ind([10 10], 1, 4)); CA1corr(:,:,sub2ind([10 10], 1, 6)); ...
    CA1corr(:,:,sub2ind([10 10], 2, 3)); CA1corr(:,:,sub2ind([10 10], 2, 5)); CA1corr(:,:,sub2ind([10 10], 4, 5))];
data = [temp(:); temp2(:)];
numel = length(temp(:));
group = [ones(numel, 1); 2*ones(numel, 1)];

[p, tbl, stats] = kruskalwallis(data, group, 'off');

fig = figure;
fig.Position = [10 10 800 800];
subplot(4, 3, 1)
hold on
histogram(temp, histBins, 'Normalization', 'probability', 'FaceColor', colors(1,:))
histogram(temp2, histBins, 'Normalization', 'probability', 'FaceColor', colors(2,:))
box off
ylim([0 0.25])
xlabel('Correlation Coefficient')
ylabel('Frequency')
title(['p = ' num2str(p)])
legend('Same Dir', 'Opposite Dir')

figure(fig)
subplot(4, 3, [4 7])
h = qqplot(temp(:));
hold on
h2 = qqplot(temp2(:));

set(h(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(1,:));
set(h(2), 'lineWidth', 2, 'color', colors(1,:))
set(h(3), 'lineWidth', 2, 'color', colors(1,:))
set(h2(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(2,:));
set(h2(2), 'lineWidth', 2, 'color', colors(2,:))
set(h2(3), 'lineWidth', 2, 'color', colors(2,:))
%% try same for CA3 same and opposite

temp = [CA3corr(:,:,sub2ind([10 10], 1, 3)); CA3corr(:,:,sub2ind([10 10], 1, 5)); CA3corr(:,:,sub2ind([10 10], 3, 5)); ...
    CA3corr(:,:,sub2ind([10 10], 2, 4)); CA3corr(:,:,sub2ind([10 10], 2, 6)); CA3corr(:,:,sub2ind([10 10], 4, 6))];

temp2 = [CA3corr(:,:,sub2ind([10 10], 1, 2)); CA3corr(:,:,sub2ind([10 10], 1, 4)); CA3corr(:,:,sub2ind([10 10], 1, 6)); ...
    CA3corr(:,:,sub2ind([10 10], 2, 3)); CA3corr(:,:,sub2ind([10 10], 2, 5)); CA3corr(:,:,sub2ind([10 10], 4, 5))];
data = [temp(:); temp2(:)];
numel = length(temp(:));
group = [ones(numel, 1); 2*ones(numel, 1)];

[p, tbl, stats] = kruskalwallis(data, group, 'off');

figure(fig)
subplot(4, 3, 2)
hold on
histogram(temp, histBins, 'Normalization', 'probability', 'FaceColor', colors(3,:))
histogram(temp2, histBins, 'Normalization', 'probability', 'FaceColor', colors(4,:))
box off
ylim([0 0.25])
xlabel('Correlation Coefficient')
ylabel('Frequency')
title(['p = ' num2str(p)])
legend('Same Dir', 'Opposite Dir')

figure(fig)
subplot(4, 3, [5 8])
h = qqplot(temp(:));
hold on
h2 = qqplot(temp2(:));

set(h(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(3,:));
set(h(2), 'lineWidth', 2, 'color', colors(3,:))
set(h(3), 'lineWidth', 2, 'color', colors(3,:))
set(h2(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(4,:));
set(h2(2), 'lineWidth', 2, 'color', colors(4,:))
set(h2(3), 'lineWidth', 2, 'color', colors(4,:))
%% try same for LS same and opposite

temp = [LScorr(:,:,sub2ind([10 10], 1, 3)); LScorr(:,:,sub2ind([10 10], 1, 5)); LScorr(:,:,sub2ind([10 10], 3, 5)); ...
    LScorr(:,:,sub2ind([10 10], 2, 4)); LScorr(:,:,sub2ind([10 10], 2, 6)); LScorr(:,:,sub2ind([10 10], 4, 6))];

temp2 = [LScorr(:,:,sub2ind([10 10], 1, 2)); LScorr(:,:,sub2ind([10 10], 1, 4)); LScorr(:,:,sub2ind([10 10], 1, 6)); ...
    LScorr(:,:,sub2ind([10 10], 2, 3)); LScorr(:,:,sub2ind([10 10], 2, 5)); LScorr(:,:,sub2ind([10 10], 4, 5))];
data = [temp(:); temp2(:)];
numel = length(temp(:));
group = [ones(numel, 1); 2*ones(numel, 1)];

[p, tbl, stats] = kruskalwallis(data, group, 'off');

figure(fig)
subplot(4, 3, 3)
hold on
histogram(temp, histBins, 'Normalization', 'probability', 'FaceColor', colors(5,:))
histogram(temp2, histBins, 'Normalization', 'probability', 'FaceColor', colors(6,:))
box off
ylim([0 0.25])
xlabel('Correlation Coefficient')
ylabel('Frequency')
title(['p = ' num2str(p)])
legend('Same Dir', 'Opposite Dir')

figure(fig)
subplot(4, 3, [6 9])
h = qqplot(temp(:));
hold on
h2 = qqplot(temp2(:));

set(h(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(5,:));
set(h(2), 'lineWidth', 2, 'color', colors(5,:))
set(h(3), 'lineWidth', 2, 'color', colors(5,:))
set(h2(1), 'marker', '.', 'markersize', 10, 'markeredgecolor', colors(6,:));
set(h2(2), 'lineWidth', 2, 'color', colors(6,:))
set(h2(3), 'lineWidth', 2, 'color', colors(6,:))


%% Get population vector correlation with distance - CA1

colors = cbrewer('qual', 'Set1', 8, 'cubic');


temp = CA1corr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA1corr(:,:,sub2ind([10 10], 3, 3));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = CA1corr(:,:,sub2ind([10 10], 5, 5));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end

% test whether they come from the same distributions
for ii = 1:size(temp, 1)
    [~, p(ii)] = kstest2(temp(:, ii), temp3(:, ii));
end

% Plot correlation distances cond 1 3 5
f = figure;
f.Position = [10 10 800 800];
subplot(4, 4, [1 ])
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA1 ')

%% CA3 Get population vector correlation with distance

temp = CA3corr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 3, 3));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = CA3corr(:,:,sub2ind([10 10], 5, 5));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end



% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, [2 ])
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('CA3 ')
xlim([0 length(binIDs)])
ylim([-0.8 1])
%% LS Get population vector correlation with distance

temp = LScorr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = LScorr(:,:,sub2ind([10 10], 3, 3));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = LScorr(:,:,sub2ind([10 10], 5, 5));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end



% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 3)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
%% Compare region 1 corr with other conds

temp = CA1corr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA1corr(:,:,sub2ind([10 10], 1, 3));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = CA1corr(:,:,sub2ind([10 10], 1, 5));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end



% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 5)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('CA1 1-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1])
%% CA3 Compare region 1 corr with other conds

temp = CA3corr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 1, 3));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = CA3corr(:,:,sub2ind([10 10], 1, 5));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end



% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 6)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('CA3 1-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1])
%% LS Compare region 1 corr with other conds

temp = LScorr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = LScorr(:,:,sub2ind([10 10], 1, 3));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = LScorr(:,:,sub2ind([10 10], 1, 5));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end

% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 7)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
title('LS 1-1 3 5')
xlim([0 length(binIDs)])
ylim([-0.8 1])
%% CA1 Compare region 1 corr with same dir and opposite dir

temp = CA1corr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA1corr(:,:,sub2ind([10 10], 1, 2));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = CA1corr(:,:,sub2ind([10 10], 1, 6));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end

% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 9)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA1 1-1 2 6')
%% CA3 Compare region 1 corr with same dir and opposite dir

temp = CA3corr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 1, 2));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = CA3corr(:,:,sub2ind([10 10], 1, 6));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end



% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 10)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA3 1-1 2 6')

%% LS Compare region 1 corr with same dir and opposite dir

temp = LScorr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = LScorr(:,:,sub2ind([10 10], 1, 2));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = LScorr(:,:,sub2ind([10 10], 1, 6));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end



% Plot correlation distances cond 1 3 5
% f = figure;
subplot(4, 4, 11)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('LS 1-1 2 6')


%% Compare regions, jump cond1
temp = CA1corr(:,:,1);
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,1);
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = LScorr(:,:,1);
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end


% f = figure;
subplot(4, 4, 13)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('compare reg.')
%% Compare regions, jump cond2
temp = CA1corr(:,:,sub2ind([10 10], 2, 2));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 2, 2));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = LScorr(:,:,sub2ind([10 10], 2, 2));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end


% f = figure;
subplot(4, 4, 14)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('compare reg. c2')


%% Compare regions, linear
temp = CA1corr(:,:,sub2ind([10 10], 9, 9));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 9, 9));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = LScorr(:,:,sub2ind([10 10], 9, 9));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end




% f = figure;
subplot(4, 4, 15)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('compare reg. c9')
%% Compare regions, linear cond7
temp = CA1corr(:,:,sub2ind([10 10], 7, 7));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 7, 7));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

temp3 = LScorr(:,:,sub2ind([10 10],  7, 7));
temp3 = triu(temp3);
temp3(~(triu(ones(size(temp3))))) = nan;
for ii = 2:length(temp3)
    temp3(ii,:) = circshift(temp3(ii,:), -ii+1);
end




% f = figure;
subplot(4, 4, 16)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp3, options)
box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('compare reg. c7')

%% CA1 Compare linear
temp = CA1corr(:,:,sub2ind([10 10], 7, 7));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA1corr(:,:,sub2ind([10 10], 9, 9));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end


f = figure;
f.Position = [10 10 800 800];
subplot(4, 4, 1)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA1 7 9')
%% CA3 Compare linear
temp = CA3corr(:,:,sub2ind([10 10], 7, 7));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 9, 9));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end


% f = figure;
subplot(4, 4, 2)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA3 7 9')

%% LS Compare linear
temp = LScorr(:,:,sub2ind([10 10], 7, 7));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = LScorr(:,:,sub2ind([10 10], 9, 9));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end

% f = figure;
subplot(4, 4, 3)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('LS 7 9')
%% CA1 Compare linear
temp = CA1corr(:,:,sub2ind([10 10], 9, 9));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA1corr(:,:,sub2ind([10 10], 9, 10));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end


% f = figure;
subplot(4, 4, 5)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA1 9-9 10')
%% CA3 Compare linear
temp = CA3corr(:,:,sub2ind([10 10], 9, 9));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = CA3corr(:,:,sub2ind([10 10], 9, 10));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end


% f = figure;
subplot(4, 4, 6)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('CA3 9-9 10')
%% LS Compare linear
temp = LScorr(:,:,sub2ind([10 10], 9, 9));
temp = triu(temp);
temp(~(triu(ones(size(temp))))) = nan;
for ii = 2:length(temp)
    temp(ii,:) = circshift(temp(ii,:), -ii+1);
end

temp2 = LScorr(:,:,sub2ind([10 10], 9, 10));
temp2 = triu(temp2);
temp2(~(triu(ones(size(temp2))))) = nan;
for ii = 2:length(temp2)
    temp2(ii,:) = circshift(temp2(ii,:), -ii+1);
end


% f = figure;
subplot(4, 4, 7)
options.handle = f;
options.error = 'std';
options.line_width = 2;
options.color_line = colors(1,:);
options.color_area = (colors(1,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
options.alpha      = 0.2;

plot_areaerrorbar(temp, options)
options.color_line = colors(2,:);
options.color_area = (colors(2,:)*1.5);
options.color_area = options.color_area/max(options.color_area);
hold on
plot_areaerrorbar(temp2, options)
options.color_line = colors(3,:);
options.color_area = (colors(3,:)*1.5);
options.color_area = options.color_area/max(options.color_area);

box off
xticks([0:7.5:40]); xticklabels([0:30:160])
xlabel('distance (cm)')
ylabel('Population Corr.')
% legend('CA1', 'CA3', 'LS')
xlim([0 length(binIDs)])
ylim([-0.8 1])
title('LS 9-9 10')

