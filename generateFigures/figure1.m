% generate position figures
clear all, close all, clc
dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

load([basename '.behavior.mat']);

colors = cbrewer('qual', 'Set2', 8, 'cubic');

%%
plotCondTrials(behavior, 'dim', "xz",  'trialType', "jump")
fig = figure(1);
fig.Position = [10, 10, 800, 800];


%% Pick one and instead plot the total position


fig = figure;
cond = 1;
subplot(5, 2, cond)
hold on
trials = find(behavior.events.trialConditions == cond);
for ii = trials
    % get timing of trials
    % then get position and plot
    trialPoints = [behavior.events.trials{ii}.x behavior.events.trials{ii}.y];
    trialPoints(end+1, :) = behavior.events.jumpLoc(ii,[1 3]);
    trialPoints(end+1, :) = behavior.events.jumpLoc(ii,[2 4]);
    [~, trialPos] = linearize(behavior, 'points', trialPoints);
    trialPos = trialPos/1000;
    scatter(trialPos(1:end-2), behavior.events.trials{ii}.z/1000, 5, 'k', 'filled')
    scatter(trialPos(1), behavior.events.trials{ii}.z(1)/1000, 10, 'g', 'filled')
    scatter(trialPos(end-2), behavior.events.trials{ii}.z(end)/1000, 10, 'r', 'filled')
    
end

xlim([- 0.1 1.8])
ylim([0.55 0.75])
axis equal
xlabel('x')
ylabel('z')
title(['ntrial = ' num2str(length(trials))])
fig.Position = [10, 10, 800, 800];

cond = 7;
subplot(5, 2, cond)
hold on
trials = find(behavior.events.trialConditions == cond);
for ii = trials
    % get timing of trials
    % then get position and plot
    trialPoints = [behavior.events.trials{ii}.x behavior.events.trials{ii}.y];
    [~, trialPos] = linearize(behavior, 'points', trialPoints);
    trialPos = trialPos/1000;
    scatter(trialPos, behavior.events.trials{ii}.z/1000, 5, 'k', 'filled')
    scatter(trialPos(1), behavior.events.trials{ii}.z(1)/1000, 10, 'g', 'filled')
    scatter(trialPos(end), behavior.events.trials{ii}.z(end)/1000, 10, 'r', 'filled')
    
end
xlim([- 0.1 1.8])
ylim([0.55 0.75])
axis equal
xlabel('x')
ylabel('z')
title(['ntrial = ' num2str(length(trials))])
fig.Position = [10, 10, 800, 800];

%% Okay good
%% Now plot the acceleration

% first smooth in x, y, and z
dt = 1/120; 
T = 120;

xPos = kalmanSmooth(behavior.position.x, dt);

figure
hold on
plot(behavior.position.x)
plot(xPos)
%% Then get x and z position
yPos = kalmanSmooth(behavior.position.y, dt);

figure
hold on
plot(behavior.position.y)
plot(yPos)

% z pos
zPos = kalmanSmooth(behavior.position.z, dt);

figure
hold on
plot(behavior.position.z)
plot(zPos)

%%
ind = length(xPos);
[~, position] = linearize(behavior, 'points', [xPos(1:ind) yPos(1:ind)]);

% then get 2 seconds to 1 second after
cond = 1;
trials = find(behavior.events.trialConditions == cond);%find(mod(behavior.events.trialConditions, 2) ==0 & behavior.events.trialConditions < 7);
trialXpos = nan(length(trials), 3*T+1);
trialXvel = nan(length(trials), 3*T);
trialXacc = nan(length(trials), 3*T-1);

count = 0;
for ii = trials
    count = count + 1;
    [~, ind] = min(abs(behavior.events.jumpTime(ii,1)-behavior.timestamps));
    px = position(ind-2*120:ind + 120)/1000;
    vel = smoothdata(diff(smoothdata(px, 'gaussian', [0 3])),'gaussian', [0 3])*120;
    acc = smoothdata(diff(vel), 'gaussian', [0 3])*120;
    
    trialXpos(count,:) = px;
    trialXvel(count, :) = vel;
    trialXacc(count, :) = acc;
   
end

fig = figure;
subplot(5, 2, 1)
imagesc(trialXpos(:, 1:end-2))
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
yticks([1 length(trials)]); yticklabels([1 length(trials)])
title('Position')
ylabel('trial #')
box off
% colorbar
subplot(5, 2, 3)
imagesc(trialXvel(:, 1:end-1))
yticks([1 length(trials)]); yticklabels([1 length(trials)])
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
title('Velocity')
ylabel('trial #')
box off
% colorbar
subplot(5, 2, 5)
imagesc((trialXacc))
yticks([1 length(trials)]); yticklabels([1 length(trials)])
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
title('Acceleration')
ylabel('trial #')
box off
% colorbar
clim([-50 50])
subplot(5, 2, 7)
hold on
yyaxis right
plot(zscore(mean(trialXpos)), 'lineWidth', 2, 'color', colors(1,:))
yyaxis left
plot(zscore(mean(trialXvel)), 'lineWidth', 2, 'color', colors(2,:))
plot(zscore(mean(trialXacc)), '-', 'lineWidth', 2, 'color', colors(3,:))
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
xlabel('Time (s)')
ylabel('z-score')
xlim([0 size(trialXacc, 2)])
fig.Position = [10 10 800 800];

%% do the same for z-position

trialZpos = nan(length(trials), 3*T+1);
trialZvel = nan(length(trials), 3*T);
trialZacc = nan(length(trials), 3*T-1);

count = 0;
for ii = trials
    count = count + 1;
    [~, ind] = min(abs(behavior.events.jumpTime(ii,2)-behavior.timestamps));
    pz = zPos(ind-2*120:ind + 120);
    vel = smoothdata(diff(smoothdata(pz, 'gaussian', [0 5])),'gaussian', [0 5]);
    acc = smoothdata(diff(vel), 'gaussian', [0 5]);
    
    trialZpos(count,:) = pz;
    trialZvel(count, :) = vel;
    trialZacc(count, :) = acc;
   
end

fig = figure;
subplot(4, 2, 1)
imagesc(trialZpos(:, 1:end-2))
yticks([1 length(trials)]); yticklabels([1 length(trials)])
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
title('Position')
box off
% colorbar
subplot(4, 2, 3)
imagesc(trialZvel(:, 1:end-1))
yticks([1 length(trials)]); yticklabels([1 length(trials)])
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
title('Velocity')
box off
% colorbar
subplot(4, 2, 5)
imagesc((trialZacc))
yticks([1 length(trials)]); yticklabels([1 length(trials)])
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
title('Acceleration')
box off
% colorbar
subplot(4,2, 7)
hold on
plot(zscore(mean(trialZpos)), 'lineWidth', 2, 'color', colors(1,:))
plot(zscore(mean(trialZvel)), 'lineWidth', 2, 'color', colors(2,:))
plot(zscore(mean(trialZacc)), 'lineWidth', 2, 'color', colors(3,:))
xticks(linspace(1, 3*T-1, 4)); xticklabels(-2:1:1)
xlabel('Time (s)')
ylabel('z-score')
xlim([0 size(trialZacc, 2)])
fig.Position = [10 10 800 800];

%% Plot overlapping jumps

fig = figure;
count = 0;
for cond = [1 3 5]
    count = count + 1;
    trials = find(behavior.events.trialConditions == cond);
    
    subplot(5, 2, cond + 2)
    hold on
    
    for ii = trials
        % get timing of trials
        % then get position and plot
        trialPoints = [behavior.events.trials{ii}.x behavior.events.trials{ii}.y];
        trialPoints(end+1, :) = behavior.events.jumpLoc(ii,[1 3]);
        trialPoints(end+1, :) = behavior.events.jumpLoc(ii,[2 4]);
        [~, trialPos] = linearize(behavior, 'points', trialPoints);
        trialPos = trialPos/1000;
        scatter(trialPos(1:end-2), behavior.events.trials{ii}.z/1000, 5,  colors(count,:), 'filled')
        scatter(trialPos(1), behavior.events.trials{ii}.z(1)/1000, 15, 'g', 'filled')
        scatter(trialPos(end-2), behavior.events.trials{ii}.z(end)/1000, 15, 'r', 'filled')
        
    end
    jumpPoints = behavior.events.jumpLoc(trials, [1 3]);
    jumpPoints = [jumpPoints; behavior.events.jumpLoc(trials, [2 4])];
   [~, jumpPos] = linearize(behavior, 'points', jumpPoints);
   jumpPos = jumpPos/1000;
   land = median(jumpPos(length(jumpPos)/2 + 1:end));
   jump = median(jumpPos(1:length(jumpPos)/2));
   scatter(jump, median(behavior.events.jumpLoc(trials, 5))/1000, 30, 'k', 'filled', '^')
   scatter(land, median(behavior.events.jumpLoc(trials, 6))/1000, 30, 'k', 'filled', 'v')
   
    xlim([- 0.1 1.8])
    ylim([0.55 0.75])
    xlabel('x')
    ylabel('z')
    title(['ntrial = ' num2str(length(trials))])
    
end

fig.Position = [10, 10, 800, 800];

%% Plot overlapping jumps
behavStart = min(behavior.events.trialIntervals(behavior.events.trialConditions< 7, 1));
behavEnd = max(behavior.events.trialIntervals(behavior.events.trialConditions< 7, 2));
fig = figure;
count = 3;
for cond = [ 7]
    count = count + 1;
    trials = find(behavior.events.trialConditions == cond);
   
        preTrials = trials(behavior.events.trialIntervals(trials, 1) < behavStart);
        postTrials = trials(behavior.events.trialIntervals(trials, 1) > behavEnd);
        
      subplot(5, 2, cond + 2)
    hold on    
    for ii = postTrials
        % get timing of trials
        % then get position and plot
        trialPoints = [behavior.events.trials{ii}.x behavior.events.trials{ii}.y];
        [~, trialPos] = linearize(behavior, 'points', trialPoints);
        trialPos = trialPos/1000;
        scatter(trialPos, behavior.events.trials{ii}.z/1000, 5,  colors(count,:), 'filled')
        scatter(trialPos(1), behavior.events.trials{ii}.z(1)/1000, 15, 'g', 'filled')
        scatter(trialPos(end), behavior.events.trials{ii}.z(end)/1000, 15, 'r', 'filled')        
    end
    xlim([- 0.1 1.8])
    ylim([0.55 0.75])
    xlabel('x')
    ylabel('z')
    title(['ntrial = ' num2str(length(postTrials))])
        
    subplot(5, 2, 1)
    hold on    
    for ii = preTrials
        % get timing of trials
        % then get position and plot
        trialPoints = [behavior.events.trials{ii}.x behavior.events.trials{ii}.y];
        [~, trialPos] = linearize(behavior, 'points', trialPoints);
        trialPos = trialPos/1000;
        scatter(trialPos, behavior.events.trials{ii}.z/1000, 5,  colors(count+1,:), 'filled')
        scatter(trialPos(1), behavior.events.trials{ii}.z(1)/1000, 15, 'g', 'filled')
        scatter(trialPos(end), behavior.events.trials{ii}.z(end)/1000, 15, 'r', 'filled')
        
    end
    xlim([- 0.1 1.8])
    ylim([0.55 0.75])
    xlabel('x')
    ylabel('z')
    title(['ntrial = ' num2str(length(preTrials))])
        
end



fig.Position = [10, 10, 800, 800];
