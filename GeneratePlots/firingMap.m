function [frMap, occuMap, xBins] = firingMap (varargin)

% makes a firing rate map of cell activity
%
%%% Inputs %%%
% 
% spikes
% 
% behavior
%
% timeInterval  Time to include in map
% 
% binSize       Width of position bin
% 
% region        Brain region (Default: hpc)
% 
% smoothMethod  Smoothing type; options: 'none', 'gaussian', 'movmean'
%
% smoothRange   Smoothing window
%
% trials        Give firing rates of individual trials, or averaged together
% 
%%% Outputs %%%
%
% frMap         firing rate map of cells
% 
% occuMap       Occupancy position map

p = inputParser;
addRequired(p, 'spikes', @isstruct);
addRequired(p, 'behavior', @isstruct); 
addParameter(p, 'timeInterval', [0 20*60], @isnumeric)
addParameter(p, 'binSize', 40, @isnumeric)
addParameter(p, 'region', 'all', @ischar)
addParameter(p, 'smoothMethod', 'gaussian', @ischar)
addParameter(p, 'smoothRange', 5, @isnumeric)
addParameter(p, 'trials', false, @islogical)

parse(p, varargin{:});

spikes = p.Results.spikes;
behavior = p.Results.behavior;
timeInterval = p.Results.timeInterval; % time interval for firing rate map
dx = p.Results.binSize; % position bin size
region = p.Results.region; % brain region
smoothMethod = p.Results.smoothMethod;
smoothRange = p.Results.smoothRange;
trials = p.Results.trials;

% hack way of dealing with this
if length(smoothMethod) == length('none')
    if sum(smoothMethod == 'none') == length('none')
        smoothRange = 1;
        smoothMethod = 'movmean';
    end
end


% params
xBins = -200:dx:2000;%0:dx:1800; % -200:dx:1900 % position bins

% Linearize position data
% angleVal = bestFitLine(behavior);
% position = (behavior.position.x + abs(min(behavior.position.x)))/cos(angleVal);
startLoc = behavior.events.startEndPos(1, 3:4);
[position, linPoints] = linearize(behavior, 'points', startLoc);
position = position-linPoints;

if length(region) == 3 && sum('all' == region) == 3
    regionVec = 1:length(spikes.times);
else
    % Get cells in region
    regionVec = zeros(1, length(spikes.times));
    for ii = 1:length(spikes.times)
        if numel(spikes.region{ii}) == length(region)
            if sum(spikes.region{ii} == region) == length(region)
                
                regionVec(ii) = ii;
            end
        end
    end
    regionVec(regionVec == 0) = [];
end

if trials == false % if want average map, not per trial map
    % Initialize firing rate map
   
    frMap = zeros(length(regionVec), length(xBins)-1);
    
    occuMap = zeros(1, length(xBins)-1);
    
    
    for interval = 1:size(timeInterval, 1)
        [~, a] = min(abs(timeInterval(interval, 1)-behavior.timestamps));
        [~, b] = min(abs(timeInterval(interval, 2)-behavior.timestamps));
        [occup, ~] = histcounts(position(a:b), xBins);
        occup = smoothdata(occup*(1/behavior.samplingRate), smoothMethod, smoothRange);
        occuMap = occuMap + occup; % total occupancy time
        count = 0;
        for ncell = regionVec % for each cell
            count = count + 1;
            [status, ~] = InIntervals(spikes.times{ncell}, timeInterval(interval,:));
            spikePos = interp1(behavior.timestamps, position, spikes.times{ncell}(status));
            [N, ~] = histcounts(spikePos, xBins);
            
            frMap(count, :) = frMap(count, :) + smoothdata(N, smoothMethod, smoothRange);
        end
    end
    frMap = frMap./occuMap; % divide by total time occupancy
%     frMap = frMap./size(timeInterval, 1); % divide by trial count (why
%     did I do this?????)

else % if want per trial map
    % Initialize firing rate map
    
    frMap = zeros(size(timeInterval, 1), length(regionVec), length(xBins)-1);
    
    occuMap = zeros(size(timeInterval, 1), length(xBins)-1);
    
    
    for interval = 1:size(timeInterval, 1)
        [~, a] = min(abs(timeInterval(interval, 1)-behavior.timestamps));
        [~, b] = min(abs(timeInterval(interval, 2)-behavior.timestamps));
        [occup, ~] = histcounts(position(a:b), xBins);
        occup = smoothdata(occup*(1/behavior.samplingRate), smoothMethod, smoothRange);
        occuMap(interval,:) = occup;
        count = 0;
        for ncell = regionVec % for each cell
            count = count + 1;
            [status, ~] = InIntervals(spikes.times{ncell}, timeInterval(interval,:));
            spikePos = interp1(behavior.timestamps, position, spikes.times{ncell}(status));
            [N, ~] = histcounts(spikePos, xBins);
            frMap(interval, count, :) = smoothdata(N, smoothMethod, smoothRange)./occup;
        end
    end
    
end

