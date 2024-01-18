function [phaseMap, xBins] = getPhaseMap(varargin)
%
% INPUTS
%  behavior
%
% spikes
%
% lfp
%
%
% OUTPUTS
%
% phaseMap
%
% varMap - variance in phase at each bin
%
% xBins
%

p = inputParser;
addRequired(p, 'behavior', @isstruct);
addRequired(p, 'spikes', @isstruct);
addRequired(p, 'lfp', @isstruct);
addRequired(p, 'Intervals', @isnumeric);
addParameter(p, 'dx', 40, @isnumeric);
addParameter(p, 'region', 'all', @ischar);

parse(p, varargin{:});

spikes = p.Results.spikes;
behavior = p.Results.behavior;
lfp = p.Results.lfp;
Intervals = p.Results.Intervals;
dx = p.Results.dx;
region = p.Results.region;

% position bins
xBins = 0:dx:1800;
xInds = 1:length(xBins)-1;

filterFreq = [6 13]/(1250*0.5);
[b, a] = butter(3, filterFreq, 'bandpass');

%
if size(Intervals, 2) > 2
    Intervals = Intervals';
end

% get list of cells
if strcmp(region, 'all')== 1
    cellIndex = 1:length(spikes.times);
else
    cellIndex = find(strcmp(spikes.region, region));
end

% get position
[position, ~, ~] = jumpPosition(behavior);

% get phase
phases = angle(hilbert(filtfilt(b, a, double(lfp.data))));

% Okay, now for each trial, get the position and phase for each cell

phaseMap = nan(length(cellIndex), length(xBins)-1);
stdMap = nan(length(cellIndex), length(xBins)-1);
%%
count = 1;
for ncell = cellIndex
    
    % get all spikes in trials
    [status, ~] = InIntervals(spikes.times{ncell}, Intervals);
    spikePos = interp1(behavior.timestamps, position, spikes.times{ncell}(status));
    spikePhase = interp1(lfp.timestamps, phases, spikes.times{ncell}(status));
    
    [~, ~,ids] = histcounts(spikePos, xBins);
    
    % go through and get average phase at each bin
    for ii = xInds
        if any(ids==ii)
            temp = spikePhase(ids==ii);
            
            phaseMap(count, ii) = circ_median(temp);
            stdMap(count, ii) = circ_std(temp);
            
        end
    end
    count = count + 1;
end

%%
% 
% figure
% subplot(2, 1, 1)
% imagesc(phaseMap)
% box off
% colorbar
% xlabel('Pos')
% ylabel('nCell')
% subplot(2, 1, 2)
% imagesc(stdMap)
% box off
% colorbar
% xlabel('Pos')
% ylabel('nCell')
