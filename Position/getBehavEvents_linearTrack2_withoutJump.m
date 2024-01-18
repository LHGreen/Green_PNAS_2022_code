% get linear session trials

clear all, close all; clc
dirName = pwd;
basename = bz_BasenameFromBasepath(dirName);

projectDir = fullfile('D:', 'BuzsakiRinzel', 'linearData');
sessionDir = fullfile(projectDir, basename);   
mkdir(sessionDir);

pos = load([basename '.pos']);

%%

if nanstd(pos(:,8)) < 10
    pos(:,8:10) = pos(:,8:10)*1000;
end
p = pos(:,[8 10]);

% Then get start/stop positions
figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(pos(:,8),pos(:,10),'.') % x-z is horizontal direction, y is height
axis([-800 1400 -1400 1200])%
% axis([600 2800 -1600 1200])%
title('start points')
[startx, starty] = ginput(2);
close

figure('units', 'normalized', 'outerposition', [0 0 1 1])
scatter(pos(:,8),pos(:,10),'.') % x-z is horizontal direction, y is height
axis([-800 1400 -1400 1200])%
% axis([600 2800 -1600 1200])
title('end points')
[endx, endy] = ginput(2);

close
%%
behavior = getBehavEvents_linearTrack2(pos, [startx starty endx endy]);


%%
plotCondTrials(behavior, 'dim', "xy")
plotCondTrials(behavior, 'dim', "x")
plotCondTrials(behavior, 'dim', "xz")

plotAllTrials(behavior, 'dim', "xzt")


%%
dataName = [basename '.behavior.mat'];
oldDir = cd(sessionDir);
save(dataName, 'behavior')

cd(oldDir)


