function equalBinMat = getSameNumBins(spikes, behavior, position, trials, binNum, region, cond)

% function to get matrix of mean firing rates from cells, with equal number
% of bins in each trial, though the time and position in each bin changes

% spikes        Buzcode format spikes file
% timepoints    Trial start time, jump time, land time, end time, format
%               [nTrials x 4]
% behavior 
%
% binNum

nBins = sum(binNum);
cumBin = cumsum(binNum);
regionVec = strcmp(region, spikes.region);
equalBinMat = zeros(sum(regionVec), nBins);

timePoints = [behavior.events.trialIntervals(trials, 1) ...
    behavior.events.jumpTime(trials, :) behavior.events.trialIntervals(trials, 2)];

startBins = nan(size(timePoints, 1), binNum(1)+1);
jumpBins = nan(size(timePoints, 1),  binNum(2) + 1);
stopBins = nan(size(timePoints, 1), binNum(3) + 1);

%%
for nTrial = 1:size(timePoints, 1)
    
    % start
    start = interp1(behavior.timestamps, position, timePoints(nTrial, 1:2));
    if diff(start) < 0
        start = start(2:-1:1);
        startBins(nTrial, :) = linspace(start(1), start(2), binNum(1) + 1);
%         startBins(nTrial, :) = startBins(nTrial, end:-1:1);
    else
        startBins(nTrial, :) = linspace(start(1), start(2), binNum(1) + 1);
    end
    
    % jump
    jump = interp1(behavior.timestamps, position, timePoints(nTrial, 2:3));
    if diff(jump) < 0
        jump = jump(2:-1:1);
        jumpBins(nTrial, :) = linspace(jump(1), jump(2), binNum(2) + 1);
%         jumpBins(nTrial, :) = jumpBins(nTrial, end:-1:1);
    else
        jumpBins(nTrial, :) = linspace(jump(1), jump(2), binNum(2) + 1);
    end
    
    
    % end
    stop = interp1(behavior.timestamps, position, timePoints(nTrial, 3:4));
    if diff(stop) < 0
        stop = stop(2:-1:1);
        stopBins(nTrial, :) = linspace(stop(1), stop(2), binNum(3) + 1);
%         stopBins(nTrial, :) = stopBins(nTrial, end:-1:1);
    else
        stopBins(nTrial, :) = linspace(stop(1), stop(2), binNum(3) + 1);
    end
    
end

count = 0;
for nCell = find(regionVec)
    count = count + 1;
    
    % start:jump
    [timeStatus, ~] = InIntervals(spikes.times{nCell}, timePoints(:, [1:2]));
    spikePos = interp1(behavior.timestamps, position, spikes.times{nCell}(timeStatus));
    
    for nTrial = 1:size(timePoints, 1)
        temp = histcounts(spikePos, startBins(nTrial, :));
        if mod(cond, 2) > 0
        equalBinMat(count, 1:binNum(1)) = equalBinMat(count, 1:binNum(1)) + temp;
        else
            equalBinMat(count, 1:binNum(1)) = equalBinMat(count, 1:binNum(1)) + temp(end:-1:1);
        end
    end
    
    % jump:land
    [timeStatus, ~] = InIntervals(spikes.times{nCell}, timePoints(:, [2:3]));
    spikePos = interp1(behavior.timestamps, position, spikes.times{nCell}(timeStatus));
    
    for nTrial = 1:size(timePoints, 1)
        temp = histcounts(spikePos, jumpBins(nTrial, :));
        if mod(cond, 2) > 1
        equalBinMat(count, (cumBin(1)+1):cumBin(2)) = equalBinMat(count, (cumBin(1)+1):cumBin(2)) + temp;
        else
            equalBinMat(count, (cumBin(1)+1):cumBin(2)) = equalBinMat(count, (cumBin(1)+1):cumBin(2)) + temp(end:-1:1);
        end
    end
    
    % land:end
    [timeStatus, ~] = InIntervals(spikes.times{nCell}, timePoints(:, [3:4]));
    spikePos = interp1(behavior.timestamps, position, spikes.times{nCell}(timeStatus));
    
    for nTrial = 1:size(timePoints, 1)
        temp = histcounts(spikePos, stopBins(nTrial, :));
        if mod(cond, 2) > 1
        equalBinMat(count, (cumBin(2)+1):cumBin(3)) = equalBinMat(count, (cumBin(2)+1):cumBin(3)) + temp;
        else
            equalBinMat(count, (cumBin(2)+1):cumBin(3)) = equalBinMat(count, (cumBin(2)+1):cumBin(3)) + temp(end:-1:1);
        end
    end
    
    % this should bin the spikes at the position for each trial, and add
    % them.
    % It'll be z-scored, so don't need to normalize by time or trials
    
end
    
    