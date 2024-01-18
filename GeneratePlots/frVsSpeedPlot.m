% Plot firing rate with speed
clear; close all

dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

load([basename '.behavior.mat']);
spikes = bz_GetSpikes('basepath', dirName);

projectDir = fullfile('D:', 'BuzsakiRinzel', 'Figures');
sessionDir = fullfile(projectDir, basename);


%%
% Get firing rates for each cell and condition

% Data matrix of trial num, velocity; NxBinNum

region = 'ca3';

dx = 50; % 5 cm, ~ 1 rat's length
nBins = 1800/dx;
posBins = 0:dx:1800;
% dV = 0.01;
% vBins = 0:dV:40;
% vVals = (vBins(1:end-1)+vBins(2:end))/2;

Vals = nan(length(behavior.events.trials), nBins);
accVals = nan(length(behavior.events.trials), nBins);

[frMap, occuMap, ~] = ...
    firingMap(spikes, behavior, 'binSize', dx, 'region', region, 'trials', true, 'timeInterval', behavior.events.trialIntervals);



% Get speed for each condition

% get linearized position and jump position
jumpLoc = nan(13, 2);
for cond = 1:6
    trials = find(behavior.events.trialConditions == cond);
    jumpLoc(cond, 1) = nanmedian(behavior.events.jumpLoc(trials, 1));
    jumpLoc(cond, 2) = nanmedian(behavior.events.jumpLoc(trials, 3));
    jumpLoc(cond+6, 1) = nanmedian(behavior.events.jumpLoc(trials, 2));
    jumpLoc(cond+6, 2) = nanmedian(behavior.events.jumpLoc(trials, 4));
end
jumpLoc(end, 1) = behavior.events.startEndPos(1, 3);
jumpLoc(end, 2) = behavior.events.startEndPos(1, 4);

% shift so that min is at 0
[position, linPoints] = linearize(behavior, 'points', jumpLoc);
landLoc = linPoints(7:12, :);
jumpLoc = linPoints(1:6, :);
position = position-linPoints(end);
landLoc = landLoc-linPoints(end);
jumpLoc = jumpLoc-linPoints(end);

%%
% get velocity at each position
for ntrial = 1:length(behavior.events.trials) % for each trial
    [~, trialPos] = linearize(behavior, 'points', [behavior.events.trials{ntrial}.x behavior.events.trials{ntrial}.y]);
    trialPos = trialPos-linPoints(end); % get trial position and velocity
    trialVel = abs(diff(trialPos))*behavior.samplingRate/1000; % get m/s
    trialAcc = diff(abs(diff(movmean(trialPos, 5))))*behavior.samplingRate.^2/1000; % get m/s^2
    trialAcc(end:end+2) = trialAcc(end);
    % Get the bin index of each trial position
    inds = discretize(trialPos(1:end-1), posBins);
    for posVal = 1:nBins % for each bin
        vels = trialVel(inds == posVal); % temporary velocity values are the velocities at that position
        Vals(ntrial, posVal) = median(vels); % save the median velocity at that position bin
        accs = trialAcc(inds==posVal);
        accVals(ntrial, posVal) = median(accs);
    end
    %     tempVal = histcounts2(trialVel, trialPos(1:end-1), vBins, posBins);
    %     [~, inds] = max(tempVal);
    Vals(ntrial, :) = smoothdata(Vals(ntrial,:), 'gaussian', 5); % smooth
    accVals(ntrial, :) = smoothdata(accVals(ntrial,:), 'gaussian', 5); 
end

% Plot velocity vs fr


%%
% plot for each trial, the median velocity and mean firing rate at each
% position
cmap = jet(nBins);
figure
hold on
for ntrial = [2 3 4 5 7 9] %length(behavior.events.trials)
    
    for ncell = 1:size(frMap, 2)
        
        scatter3(Vals(ntrial, :), frMap(ntrial, ncell,:), 1:nBins, 8, cmap, 'filled')
    end
    xlim([0 3])
    %     plot(linspace(0,3, 80), Vals(ntrial,:)*10)
end

xlabel('velocity (m/s)')
ylabel('Firing Rate (Hz)')
zlabel('position')

%%
figname = [basename '_region' region '_VelocityVsFr.png'];
oldDir = cd(sessionDir);
saveas(gcf, figname)
cd(oldDir)
close


% colors = zeros(8, 3);
% colors(:,3) = linspace(1, 256, 8);
%
% colors = [30, 86, 49;...
%     164, 222, 2;...
%     118, 186, 27;...
%     76, 154, 42;...
%     104, 187, 89;...
%     172, 223, 135;
%     0.5 0.5 0.5;...
%     0 0 0];

colors = [255, 244, 191;...
    88, 204, 237;...
    254, 190, 3;...
    56, 149, 211;...
    255, 126, 0;...
    7, 47, 95;...
    0.8 0.8 0.8;...
    0 0 0];

colors = colors/256;
figure('units', 'normalized', 'outerposition', [0 0 1 1])
numCells = size(frMap, 2);
for ncell = 1:numCells
    %     subplot(4,4,ncell)
    %     subplot(ceil(sqrt(numCells)), ceil(sqrt(numCells)), ncell)
    hold on
    for cond = 1
        trials = find(behavior.events.trialConditions == cond);
        for ntrial = trials
            scatter(Vals(ntrial, :), frMap(ntrial, ncell,:), 5, ...
                'MarkerFaceColor', colors(cond, :),'MarkerEdgeColor', colors(cond, :))
        end
    end
    xlim([0 3])
    ylim([0 80])
    xlabel('velocity (m/s)')
    ylabel('Firing Rate (Hz)')
end


%

figname = [basename '_region' region  '_AllCellsVelocityVsFrOB.png'];
oldDir = cd(sessionDir);
saveas(gcf, figname)
cd(oldDir)
close

% Plot per cell?
% make legend
% colors = [255, 244, 191;...
%     88, 204, 237;...
%     254, 190, 3;...
%     56, 149, 211;...
%     255, 126, 0;...
%     7, 47, 95;...
%     0.8 0.8 0.8;...
%     0 0 0];
%
% colors = colors/256;
% figure;
% hold on
% for cond = 1:6
%     scatter(cond, ones(1,1), 20, 'MarkerFaceColor', colors(cond, :),...
%         'MarkerEdgeColor', colors(cond, :));
% end
% xlabel('condition')
%
% figname = ['colorLegend.png'];
% oldDir = cd(sessionDir);
% saveas(gcf, figname)
% cd(oldDir)
% close

%
% show velocity of each trial
[~, condOrder] = sort(behavior.events.trialConditions);
figure
imagesc(Vals(condOrder, :))
xticks(1:5:size(Vals,2)+1); xticklabels(posBins(1:5:end)/10)
yticks([1 size(Vals, 1)]); yticklabels([1 length(behavior.events.trials)]);
colorbar
box off
xlabel('position (cm)')
ylabel('trial #')
figname = [basename '_velocity.png'];
oldDir = cd(sessionDir);
saveas(gcf, figname)
cd(oldDir)
close

% show mean frMap for each condition
for cond = 1:8
    ntrials = find(behavior.events.trialConditions == cond);
    figure
    imagesc(squeeze(mean(frMap(ntrials, :,:))))
    xticks(1:5:size(frMap,3)+1); xticklabels(posBins(1:5:end)/10)
    yticks([1 size(frMap, 2)]); yticklabels([1 size(frMap, 2)]);
    colorbar
    box off
    xlabel('position (cm)')
    ylabel('cell #')
    if cond < 7
        hold on
        plot((jumpLoc(cond)/dx)*ones(1, 2), [0 size(frMap, 1)], 'r', 'lineWidth', 1)
        plot((landLoc(cond)/dx)*ones(1, 2), [0 size(frMap, 1)], 'g', 'lineWidth', 1)
    end
    figname = [basename '_region_' region '_cond' num2str(cond) '_firingMap.png'];
    oldDir = cd(sessionDir);
    saveas(gcf, figname)
    cd(oldDir)
    close
end

%% plot mean velocity of each condition

meanVel = nan(8, size(Vals, 2));
semVel = nan(8, size(Vals, 2));
for ii = 1:8
    trials = find(behavior.events.trialConditions == ii);
    meanVel(ii,:) = nanmean(Vals(trials, :));
    semVel(ii,:) = sqrt(nanvar(Vals(trials, :))); %./sqrt(length(trials));
end

xx = (posBins(1:end-1)+posBins(2:end))/2000;
fig = figure;
fig.Position = [10 10 800 800];
for ii = 1:8
    subplot(4, 2, ii)
    hold on
    id = ~isnan(meanVel(ii,:));
    fill([xx(id)'; flipud(xx(id)')], [meanVel(ii,id)'-semVel(ii,id)'; flipud(meanVel(ii,id)'+semVel(ii,id)')], 'k', 'facealpha', 0.2, 'linestyle', 'none', 'HandleVisibility', 'off')
    line(xx, meanVel(ii,:), 'color', 'k', 'lineWidth', 1)
    %     plot(posBins(1:end-1)/1000, meanVel(ii,:), 'k', 'lineWidth', 1);
    if ii < 7
        %         hold on
        plot((jumpLoc(ii)/(1000))*ones(1, 2), [0 2], 'r', 'lineWidth', 2)
        plot((landLoc(ii)/(1000))*ones(1, 2), [0 2], 'g', 'lineWidth', 2)
    end
    ylim([0 2])
    box off
end
han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Velocity (m/s)');
xlabel(han,'Position (m)');


%% plot mean acceleration of each condition

meanAcc = nan(8, size(accVals, 2));
semAcc = nan(8, size(accVals, 2));
for ii = 1:8
    trials = find(behavior.events.trialConditions == ii);
    meanAcc(ii,:) = nanmean(accVals(trials, :));
    semAcc(ii,:) = sqrt(nanvar(accVals(trials, :))); %./sqrt(length(trials));
end

xx = (posBins(1:end-1)+posBins(2:end))/2000;
fig = figure;
fig.Position = [10 10 800 800];
for ii = 1:8
    subplot(4, 2, ii)
    hold on
    id = ~isnan(meanAcc(ii,:));
    fill([xx(id)'; flipud(xx(id)')], [meanAcc(ii,id)'-semAcc(ii,id)'; flipud(meanAcc(ii,id)'+semAcc(ii,id)')], 'k', 'facealpha', 0.2, 'linestyle', 'none', 'HandleVisibility', 'off')
    line(xx, meanAcc(ii,:), 'color', 'k', 'lineWidth', 1)
    %     plot(posBins(1:end-1)/1000, meanVel(ii,:), 'k', 'lineWidth', 1);
    if ii < 7
        %         hold on
        plot((jumpLoc(ii)/(1000))*ones(1, 2), [-10 10], 'r', 'lineWidth', 2)
        plot((landLoc(ii)/(1000))*ones(1, 2), [-10 10], 'g', 'lineWidth', 2)
    end
    ylim([-10 10])
    box off
end
han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Acceleration (m/s^2)');
xlabel(han,'Position (m)');


%% Plot mean velocity against mean firing rate for each position

% get indices of jump and land positions
% jumpBin = round(interp1(xx, 1:length(xx), jumpLoc/1000));


fig = figure;
for ii = 1:8
    trials = find(behavior.events.trialConditions == ii);
    subplot(4, 2, ii)
    hold on
    id = find(~isnan(meanVel(ii,:)));
    jumpBin = round(interp1(xx(id), 1:length(id), jumpLoc/1000));
    landBin = round(interp1(xx(id), 1:length(id), landLoc/1000));
    
    fr = squeeze(mean(mean(frMap(trials, :, id)), 2));
    plot(meanVel(ii, id), fr, 'k', 'lineWidth', 2)
    if mod(ii, 2)==1
    scatter(meanVel(ii, id(1)), fr(1), 10, 'r', 'filled')
    scatter(meanVel(ii, id(end)), fr(end), 10, 'g', 'filled')
    else
        scatter(meanVel(ii, id(1)), fr(1), 10, 'g', 'filled')
        scatter(meanVel(ii, id(end)), fr(end), 10, 'r', 'filled')
    end
    if ii < 7
        scatter(meanVel(ii, id(jumpBin(ii))), fr(jumpBin(ii)), 10, 'c', 'filled')
        scatter(meanVel(ii, id(landBin(ii))), fr(landBin(ii)), 10, 'm', 'filled')
    end
    xlim([0 2])
    ylim([0.5 2.2])
    box off
end

han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Velocity (m/s)');
ylabel(han,'Firing Rate (Hz)');


%% Plot mean acceleration against mean firing rate for each position

% get indices of jump and land positions
% jumpBin = round(interp1(xx, 1:length(xx), jumpLoc/1000));


fig = figure;
for ii = 1:8
    trials = find(behavior.events.trialConditions == ii);
    subplot(4, 2, ii)
    hold on
    id = find(~isnan(meanAcc(ii,:)));
    jumpBin = round(interp1(xx(id), 1:length(id), jumpLoc/1000));
    landBin = round(interp1(xx(id), 1:length(id), landLoc/1000));
    
    fr = squeeze(mean(mean(frMap(trials, :, id)), 2));
    plot(meanAcc(ii, id), fr, 'k', 'lineWidth', 2)
    if mod(ii, 2)==1
    scatter(meanAcc(ii, id(1)), fr(1), 10, 'r', 'filled')
    scatter(meanAcc(ii, id(end)), fr(end), 10, 'g', 'filled')
    else
        scatter(meanAcc(ii, id(1)), fr(1), 10, 'g', 'filled')
        scatter(meanAcc(ii, id(end)), fr(end), 10, 'r', 'filled')
    end
    if ii < 7
        scatter(meanAcc(ii, id(jumpBin(ii))), fr(jumpBin(ii)), 10, 'c', 'filled')
        scatter(meanAcc(ii, id(landBin(ii))), fr(landBin(ii)), 10, 'm', 'filled')
    end
    xlim([-8 8])
    ylim([0.5 3])
    box off
end

han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Acceleration (m/s^2)');
ylabel(han,'Firing Rate (Hz)');


