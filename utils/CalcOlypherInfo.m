function [trackInfo] = CalcOlypherInfo(firingRates, posIDs, numBins)

% USAGE
% [track_info] = olypherInfo(data)
% 
% INPUTS
%     firingRates - Ncells, x T matrix of all spike times for all cells
%     
%     posIDs - T vector of bin of each position
%
%     xBins - binCenters
%     
%     
% OUTPUTS
%     track_info - Ncells x nBins matrix of information scores across all cells and bins
%     
% TO-DO
% Add varargin stuff
%
% Add making the firing rates and position data part of the code
%       Laura Green, 2021
%% Info Analysis

numCells = size(firingRates, 1);
numTime = size(firingRates, 2);

trackInfo = zeros(numCells, numBins);

for nCell = 1:numCells % for all cells
    for tt = 1:numTime % for all time bins
            
            k = firingRates(nCell, tt);
        
%             k = frMap(nCell, nBin, nTrial); % number of spikes at bin and trial number
            temp = posIDs==posIDs(tt); % all time intervals with position x
            pKx = (sum(firingRates(nCell, temp) == k))/sum(temp); % probability of observing response k at location x
            pK = (sum(firingRates(nCell, :) == k))/numTime; % probability of observing response k
            
            if pK == 0 || pKx == 0 || pKx < pK
               trackInfo(nCell, posIDs(tt)) = trackInfo(nCell, posIDs(tt));
                
            else
%                 fprintf(['ncell = ' num2str(nCell) ', bin = ' num2str(posIDs(tt)) ', tt = ' num2str(tt) '\n'])
                
                trackInfo(nCell, posIDs(tt)) = trackInfo(nCell, posIDs(tt)) + (pKx*log2(pKx./pK));
%                 fprintf([num2str(pKx*log2(pKx./pK)) '\n'])
            end
     
    end
end



