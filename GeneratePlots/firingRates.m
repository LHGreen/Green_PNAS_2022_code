function [frRates, timeBins] = firingRates(varargin)
% function to get firing rates in time
%
% Usage
%
% [frRates, timeBins] = firingRates(varargin)
%
% Inputs
%
% spikes        buzcode format .spikes struct
%
% timeInterval  Interval of spikes to use
%
% refTime       Reference time for 0 s of spikes
%
% Outputs
%
% firing Rates  Nxnt size matrix, N = number of cells
% 
% timeBins      bins used to get firing rates

% Should probably edit this so that the interval is given, not just
% automatically -1.5 to 1 s

p = inputParser;
addRequired(p, 'spikes', @isstruct);
% addRequired(p, 'behavior', @isstruct);
addRequired(p, 'timeInterval', @isnumeric);
addRequired(p, 'refTime', @isnumeric);
addParameter(p, 'dt', 0.05, @isnumeric)
addParameter(p, 'region', 'all', @ischar)
addParameter(p, 'trials', false, @islogical)
addParameter(p, 'winSize', 10, @isnumeric)

parse(p, varargin{:});

spikes = p.Results.spikes;
% behavior = p.Results.behavior;
timeInterval = p.Results.timeInterval;
refTime = p.Results.refTime;
dt = p.Results.dt;
region = p.Results.region;
trials = p.Results.trials;
winSize = p.Results.winSize;

%
timeBins = -1.5:dt:1; % time bins 
if strcmp(region, 'all') 
    regionVec = 1:length(spikes.times);      
else
    regionVec = find(strcmp(spikes.region, region));
end

% Get firing rates of the neurons
if trials
    frRates = zeros(length(regionVec), length(timeBins)-1, size(timeInterval, 1));
else
    frRates = zeros(length(regionVec), length(timeBins)-1);
end

% if ~isempty(regionVec) 
if trials
    for trial = 1:size(timeInterval, 1) % for each time interval
        count = 0;
        for neuron = regionVec % for each neuron
            count = count + 1;
            times = spikes.times{neuron} - refTime(trial);
            [N, ~] = histcounts(times, timeBins); % Get count of spikes in each time bin
            N = smoothdata(N, 'gaussian', winSize);
            frRates(count,:, trial) = frRates(count,:, trial) + N;
            
        end
    end
    
    frRates = frRates./dt; % convert to spikes/second
    
else
    for trial = 1:size(timeInterval, 1) % for each time interval
        count = 0;
        for neuron = regionVec % for each neuron
            count = count + 1;
            times = spikes.times{neuron} - refTime(trial);
            [N, ~] = histcounts(times, timeBins); % Get count of spikes in each time bin
            N = smoothdata(N, 'gaussian', winSize);
            frRates(count,:) = frRates(count,:) + N;
            
        end
    end
    frRates = frRates/size(timeInterval, 1); % average across trials
    frRates = frRates./dt; % convert to spikes/second
    
end

        
    