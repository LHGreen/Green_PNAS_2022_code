% Make presentation figures
dx = 10; % in mm
posBins = -400:dx:820;
nbins = length(posBins)-1;
trackLength = 1220;

% Save info
projectDir = fullfile('H:', 'Personal', 'BuzsakiRinzel', 'Figures');
dirName = pwd;

% get ratemaps from the remapping code

% to make
% - firing rates for linear and jump trials
%% Get indices of silent cells
ActiveCells = find(mean(rateMapLinear(:,:,1)') ~= 0);

eCells = find(mean(rateMapLinear(:,:,1)') < 20);

[~,~, Order] = sort_cells(zscore(rateMapLinear(:, 1:120, 1)')');

ActiveOrder = intersect(Order, ActiveCells, 'stable');

[~,~, jumpOrder] = sort_cells(zscore(rateMapJump(:,1:120, 1)')');
%% Firing rates

% linear

for ii = 1:2
%     subplot(1,2,ii)
    figure
    imagesc(sort_cells(zscore(rateMapLinear(ActiveCells, 1:120, ii)')'))
    caxis([-2.5 2.5])
    xticks(1:30:120); xticklabels(round((1:30:120)/cos(angleVal)))
    xlabel('Position (cm)')
    ylabel('Cell')
    h = colorbar;
    ylabel(h, 'z-score firing rate')
    title(['Direction ' num2str(ii)])
    figname = ['CA3rateMaplinearDir' num2str(ii) '.png'];
    cd(projectDir)
    saveas(gcf, figname)
    cd(dirName)

end


% figure
for ii = 1:6
%     subplot(3,2,ii)
figure
    imagesc(sort_cells(zscore(rateMapJump(ActiveCells, 1:120, ii)')'))
    caxis([-2.5 2.5])
    xticks(1:30:nbins); xticklabels(round((1:30:120)/cos(angleVal)))
    xlabel('Position (cm)')
    ylabel('Cell')
    title(['Condition ' num2str(ii)])
    h = colorbar;
    ylabel(h, 'z-score firing rate')
    figname = ['CA3rateMapJumpCond' num2str(ii) '.png'];
    cd(projectDir)
    saveas(gcf, figname)
    cd(dirName)
%     close
end


% Remapping
% figure
for ii = 1:2
%     subplot(1,2,ii)
figure
    imagesc((zscore(rateMapLinear(ActiveOrder, 1:120, ii)')'))
    caxis([-2.5 2.5])
    xticks(1:30:nbins); xticklabels(round((1:30:120)/cos(angleVal)))
    xlabel('Position (cm)')
    ylabel('Cell')
    h = colorbar;
    ylabel(h, 'z-score firing rate')
    title(['Direction ' num2str(ii)])
    figname = ['CA3OrderLinearDir' num2str(ii) '.png'];
    cd(projectDir)
    saveas(gcf, figname)
    cd(dirName)
end

% figure
for ii = 1:6
%     subplot(3,2,ii)
figure
    imagesc((zscore(rateMapJump(ActiveOrder, 1:120, ii)')'))
    caxis([-2.5 2.5])
    xticks(1:30:nbins); xticklabels(round((1:30:120)/cos(angleVal)))
    xlabel('Position (cm)')
    ylabel('Cell')
    h = colorbar;
    ylabel(h, 'z-score firing rate')
    title(['Condition ' num2str(ii)])
    figname = ['CA3OrderJumpCond' num2str(ii) '.png'];
    cd(projectDir)
    saveas(gcf, figname)
    cd(dirName)
end

% Remapping
% figure
for ii = 1:2
%     subplot(1,2,ii)
figure
    imagesc((zscore(rateMapLinear(jumpOrder, 1:120, ii)')'))
    caxis([-2.5 2.5])
    xticks(1:30:nbins); xticklabels(round((1:30:120)/cos(angleVal)))
    xlabel('Position (cm)')
    ylabel('Cell')
    h = colorbar;
    ylabel(h, 'z-score firing rate')
    title(['Direction ' num2str(ii)])
    figname = ['CA3OrderCond1LinearDir' num2str(ii) '.png'];
    cd(projectDir)
    saveas(gcf, figname)
    cd(dirName)
end

% linear
% figure
for ii = 1:6
%     subplot(3,2,ii)
figure
    imagesc((zscore(rateMapJump(jumpOrder, 1:120, ii)')'))
    caxis([-2.5 2.5])
    xticks(1:30:nbins); xticklabels(round((1:30:120)/cos(angleVal)))
    xlabel('Position (cm)')
    ylabel('Cell')
    h = colorbar;
    ylabel(h, 'z-score firing rate')
    title(['Condition ' num2str(ii)])
    figname = ['CA3OrderCond1JumpCond' num2str(ii) '.png'];
    cd(projectDir)
    saveas(gcf, figname)
    cd(dirName)
end

%% Get trials in dir 1

x = nan(sum(behavior.events.trialConditions == 1), 200);
y = nan(sum(behavior.events.trialConditions == 1), 200);
count = 1;
for ii = find(behavior.events.trialConditions == 1)
    x(count,:) = behavior.events.trials{ii}.x;
    y(count,:) = behavior.events.trials{ii}.y;
    count = count+1;
end

%%
figure; plot(x', y', 'k', 'lineWidth', 1)
hold on;
scatter(x(:,1), y(:,1), 5,'g', 'filled')
scatter(x(:,end), y(:,end), 5,'r', 'filled')
xticks([])
yticks([]);


%% 

figure
scatter(cellPosxJump{1,1}/cos(cosAngle), 1:length(cellPosxJump{1,1}))

%%
figure
plot(movmean(double(lfp.data(timeStatus{2})), 7), 'k', 'lineWidth', 1)
xlabel('time')
ylabel('amplitude')
set(gca,'TickDir','out'); 
ylim([-4000 4000])
yticklabels('')
xticks(125:1250:sum(timeStatus{2})); xticklabels(-1:1)
axis off

%%

[status, ~,~] = InIntervals(lfp.timestamps, behavior.events.trialIntervals(5,:));
figure
plot(movmean(double(lfp.data(status)), 7), 'k', 'lineWidth', 1)
xlabel('time')
ylabel('amplitude')
set(gca,'TickDir','out'); 
ylim([-4000 4000])
yticklabels('')
xticks(125:1250:sum(timeStatus{2})); xticklabels(-1:1)
axis off

%% Plot firing rate of that interneurons

figure
plot(allRateMapJump{6}(26,:,6), 'k', 'lineWidth', 1)
xticks(0:20:120); xticklabels(round((0:20:120)/cos(0.81), -1))
ylabel('Firing Rate (Hz)')
xlabel('Position (cm)')
set(gca,'TickDir','out'); 
box off


