% put together the linear track trials and jump trials
clear all, close all; clc
dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

% load([basename '.jumpPiezo.behavior.mat']);   
% 
% linear = load([basename '.linearTrack.behavior.mat']);

projectDir = fullfile('D:', 'BuzsakiRinzel', 'Figures');
sessionDir = fullfile(projectDir, basename);   
mkdir(sessionDir);

saveDir = fullfile('D:', 'BuzsakiRinzel', 'Data', basename);

pos = load([basename '.pos']);

%%
behavior = bz_getJumpBehav_NoPiezoEditVersion(pos);
% behavior = bz_getJumpBehav_piezoV2(pos, 2);
% behavior = bz_getJumpBehav_piezo(pos, 2);

plotCondTrials(behavior, 'dim', "xy", 'trialType', "jump")
plotCondTrials(behavior, 'dim', "x")
plotCondTrials(behavior, 'dim', "xz",  'trialType', "jump")

%%
figure
scatter(diff(behavior.events.jumpTime'), behavior.events.trialConditions, 20, 'b', 'filled')
xlabel('air time')
ylabel('trial condition')

% histogram of jump/land locations
fig = figure;
for ii = 1:6
    try
    trials = find(behavior.events.trialConditions == ii);
    subplot(6, 2, ii*2-1)
    bins = mean(behavior.events.jumpLoc(trials, 1))-80:10:mean(behavior.events.jumpLoc(trials, 1))+80;
    histogram(behavior.events.jumpLoc(trials, 1), bins)
    box off
    title(['cond ' num2str(ii) ' jump'])
    subplot(6, 2, ii*2)
    bins = mean(behavior.events.jumpLoc(trials, 2))-80:10:mean(behavior.events.jumpLoc(trials, 2))+80;
    histogram(behavior.events.jumpLoc(trials, 2), bins)
    title(['cond ' num2str(ii) ' land'])
    box off
    catch
        
    end
end
fig.Position = [10 10 800 800];

%% get rid of out of range trials
% indexVal = [];
% for ii = 1:length(behavior.events.trials)
%     if mod(behavior.events.trialConditions(ii), 2) == 1
%         if ~isempty(find(behavior.events.trials{ii}.x < behavior.events.startEndPos(1, 1)-5, 1))
%             indexVal = [indexVal ii];
%             
%         end
%     end
% end
%%
indexVal = [33];
for ii = indexVal
    behavior.events.trials(ii) = [];
    behavior.events.trialConditions(ii)= [];
    behavior.events.trialIntervals(ii,:) = [];
    behavior.events.jumpTime(ii,:) = [];
    behavior.events.jumpLoc(ii,:) = [];
end
%%

linear = getBehavEvents_linearTrack2(pos, behavior.events.startEndPos);

plotCondTrials(linear, 'dim', "xy")
plotCondTrials(linear, 'dim', "x")
plotCondTrials(linear, 'dim', "xz")

plotAllTrials(linear, 'dim', "xzt")


%%
for ii = 1:length(linear.events.trials)
    behavior.events.trials{end+1} = linear.events.trials{ii};
end

behavior.events.trialConditions = [behavior.events.trialConditions linear.events.trialConditions+6];

behavior.events.trialIntervals = [behavior.events.trialIntervals; linear.events.trialIntervals];

behavior.events.jumpTime = [behavior.events.jumpTime; nan(length(linear.events.trials), 2)];

behavior.events.jumpLoc = [behavior.events.jumpLoc; nan(length(linear.events.trials), 6)];

behavior.events.conditionType{end+1} = 'linearTrack';

fprintf('done \n')
%%
plotCondTrials(behavior, 'dim', "xy")
plotCondTrials(behavior, 'dim', "x")
plotCondTrials(behavior, 'dim', "xz")


%% save

% cd(saveDir)

filename = [basename '.behavior.mat'];
save(filename, 'behavior')
% cd(dirName)

%% Get jumploc in yz if necessary
% behavior.events.jumpLoc = [behavior.events.jumpLoc(:,1:2) behavior.events.jumpLoc];
% 
% trials = find(behavior.events.trialConditions < 7);
% for ii = trials
%     [~, a] = min(abs(behavior.events.trials{ii}.timestamps-behavior.events.jumpTime(ii,1)));
%     [~, b] = min(abs(behavior.events.trials{ii}.timestamps-behavior.events.jumpTime(ii,2)));
%     behavior.events.jumpLoc(ii, 3:4) = [behavior.events.trials{ii}.y(a) behavior.events.trials{ii}.y(b)];
%     behavior.events.jumpLoc(ii, 5:6) = [behavior.events.trials{ii}.z(a) behavior.events.trials{ii}.z(b)];
% end
    
%%
% histogram of jump/land locations
figure
for ii = 1:6
    try
    trials = find(behavior.events.trialConditions == ii);
    subplot(6, 2, ii*2-1)
    bins = mean(behavior.events.jumpLoc(trials, 1))-60:10:mean(behavior.events.jumpLoc(trials, 1))+60;
    histogram(behavior.events.jumpLoc(trials, 1), bins)
    box off
    title(['cond ' num2str(ii) ' jump'])
    subplot(6, 2, ii*2)
    bins = mean(behavior.events.jumpLoc(trials, 2))-60:10:mean(behavior.events.jumpLoc(trials, 2))+60;
    histogram(behavior.events.jumpLoc(trials, 2), bins)
    title(['cond ' num2str(ii) ' land'])
    box off
    catch
        
    end
end

sgtitle('jump-land location (mm)')
figname = 'jumpLocationHistogram.png';
oldDir = cd(sessionDir);
saveas(gcf, figname)
cd(oldDir)
close                                                                                                                       

% x z plots for all conditions
plotCondTrials(behavior, 'dim', "xz",  'trialType', "jump")
figname = 'xzPos.png';
oldDir = cd(sessionDir);
saveas(gcf, figname)
cd(oldDir)
close

% x y plots for all conditions
plotCondTrials(behavior, 'dim', "xy",  'trialType', "jump")
figname = 'xyPos.png';
oldDir = cd(sessionDir);
saveas(gcf, figname)
cd(oldDir)
close

