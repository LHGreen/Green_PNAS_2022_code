function [infoRate, infoRateAtPos] = getInfoRate(varargin)

%%%
%  Output
%  
%  infoRate:      Information rate of each cell
%  
%  infoRateAtPos:     Information rate of each cell at each position
%  
%  Input
%  
%  behavior
%  
%  spikes
%  
%  binSize
%
%  trials (?)
 
 %% get inputs
 
 p = inputParser();
 addRequired(p, 'behavior', @isstruct);
 addRequired(p, 'spikes', @isstruct);
 addParameter(p, 'binSize', 2.5, @isnumeric);
 addParameter(p, 'trials', [], @isnumeric);
 
 parse(p, varargin{:})
 behavior = p.Results.behavior;
 spikes = p.Results.spikes;
 binSize = p.Results.binSize;
 trials = p.Results.trials;
 
 %% Make bins
 bins = 0:binSize:180;
 
 
 %% get mean firing rates
 ncells = length(spikes.times);
 meanFr = nan(ncells, 1);
 
 for ii = 1:ncells
     meanFr(ii) = length(spikes.times{ii})/(behavior.timestamps(end)); % total # spikes/total time of recording
 end
 
 %% Get firing maps and occupancy
 
 [frMap, occuMap, xBins] = firingMap(spikes, behavior, ...
     'timeInterval', behavior.events.trialIntervals(behavior.events.trialConditions < 7,:), 'binSize', binSize);
 
 p_i = occuMap./sum(occuMap); % probability of occupying that position
 
 nBins = length(xBins)-1;
 infoRateAtPos = nan(ncells, nBins);
 
 
 frMapRatio = log2(frMap./meanFr);
 frMapRatio(isinf(frMapRatio)) = 0;
 infoRateAtPos = p_i.*frMap.*frMapRatio;
 
 infoRate = nansum(infoRateAtPos, 2);
 
 