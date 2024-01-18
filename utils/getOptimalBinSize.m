function [binSize, cost] = getOptimalBinSize(data, T, ntrials, binRange, plotCost, findMin)
%
% Usage
% [binSize, cost] = getOptimalBinSize(data, T, ntrials, binRange, plotCost, findMin)
%
% INPUTS
%
% data - any time series data, organized into [timesteps x 1] format
% (concatenate trials)
%
% ntrials - number of trials in data
%
% T - duration of data
%
% binRange - range of bin widths to try
%
% plotCost - option to plot cost function (true or false)
%
% findMin - find min of cost function or get elbow of the curve(true or
% false)
%
% OUTPUTS
%
% binSize - optimal bin size for creating histograms or PSTHs, according
% to algorithm in Shimazaki and Shinomoto (2007)
% 
% cost - cost function associated with each bin size - you may want to
% smooth the cost function and choose a different bin size

%%
N = ceil(T./binRange); % number of bins in the time series

n = ntrials; % number of trials

costFunction = @(mean_k, var_k, binSize, n) (2*mean_k-var_k)./(n*binSize).^2;

cost = nan(length(binRange), 1);

for ii = 1:length(binRange)
    psth = histcounts(data(:), N(ii));
    
    mean_k = mean(psth);
    var_k = var(psth);
    
    cost(ii) = costFunction(mean_k, var_k, binRange(ii), n);
    
end



if findMin
    [~, bestBinID] =min(cost);
else % find elbow of curve (code from Jonah on stackoverflow)
    % https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
    % There is a trade-off here, for most purposes, it seems fine.
    nPoints = length(cost);
    allCoord = [1:nPoints; cost']';
    
    firstPoint = allCoord(1,:);
    lineVec = allCoord(end, :)-firstPoint;
    lineVecN = lineVec/sqrt(sum(lineVec.^2));
    
    vecFromFirst = bsxfun(@minus, allCoord, firstPoint);
    
    scalarProduct = dot(vecFromFirst, repmat(lineVecN, nPoints, 1), 2);
    
    vecFromFirstParallel = scalarProduct * lineVecN;
    
    vecToLine = vecFromFirst - vecFromFirstParallel;
    
    distToLine = sqrt(sum(vecToLine.^2, 2));
    
    [~, bestBinID] = max(distToLine);
end

binSize = binRange(bestBinID);

if plotCost
    figure
    plot(binRange, cost, 'k', 'lineWidth', 2)
    hold on
    scatter(binSize, cost(bestBinID), 30, 'r', 'filled')
    box off
    xlabel('Bin Size')
    ylabel('Cost')
    title(['Optimal Bin Size = ' num2str(binSize)])
    
end

%% run some test code here to get spike times organized by trials
% 
% trials = find(behavior.events.trialConditions ==2);
% data = [];
% for ntrial = 1:length(trials)
%     [status, ~, ~] = InIntervals(spikes.times{26}, behavior.events.trialIntervals(trials(ntrial), :));
%     
%     data = [data; spikes.times{26}(status)-behavior.events.trialIntervals(trials(ntrial), 1)];
% end
    