function [bestWidth] = getBestSmoothingWidth(data, widths, nTrials, interval, plotCost, findMin)
% 
% USAGE
% [bestWidth] = getBestSmoothingWidth(data, widths, interval)
% 
% INPUTS
% 
% data - nSpikes x 1 vector of all spikeTimes
% 
% widths = range of bin widths to tests
%
% interval = 2x1 vector of start and end times
%
% ntrials = number of trials
% 
% OUTPUTS
% 
% bestWidth - optimal smoothing width
%
% From Cunningham 2011 paper on optimal bin size
%% run cost function

a = interval(1);
b = interval(2);
n = nTrials;

psi = @(t, w, a, b) (1/(sqrt(pi)*4*w)) * exp(-(t^2./(4*w.^2))) * (erf((2*b + t)./(2*w))-erf((2*a+t)./(2*w)));
% t = time
% w = bin width
% a = start time
% b = end time


cost = nan(length(widths), 1);
kernelMat = nan(length(data));
psiMat = nan(length(data));

for ww = 1:length(widths)
    kernel = gausswin(widths(ww));
    for ii = 1:length(data)
        for jj = ii + 1:length(data)
            psiMat(ii, jj) = psi(data(ii)-data(jj), widths(ww), a, b);
            
            kernelMat(ii, jj) = kernel*(data(ii)-data(jj));
        end
    end
    
    cost(ww) = -(4/n^2)*sum(triu(kernelMat, 1)) + (1/n^2)*sum(psi);
end

if findMin
 [~, bestWidthID] = min(cost);
else
% find elbow of curve (code from Jonah on stackoverflow)
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
    
    [~, bestWidthID] = max(distToLine);
end

bestWidth = widths(bestWidthID);

if plotCost
    figure
    plot(widths, cost, 'k', 'lineWidth', 2)
    hold on
    scatter(bestWidth, cost(bestWidthID), 30, 'r', 'filled')
    box off
    xlabel('kernel width')
    ylabel('Cost')
    title(['Optimal Bin Size = ' num2str(bestWidth)])
    
end


