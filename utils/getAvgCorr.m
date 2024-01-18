function [mat] = getAvgCorr(A, sym)
% get average correlation data matrix

nBins = length(A);


if sym
    temp = nan(nBins, nBins);
    
    for ii = 1:nBins
        vals = diag(A, ii-1); % get the diagonal
        temp(1:length(vals), ii) = vals;
    end
    
else
    
    temp = nan(nBins*2, nBins);
    
    for ii = 1:nBins
        vals = [diag(A, ii-1); diag(A, -(ii-1))]; % get the diagonal and all off-diagonals
        temp(1:length(vals), ii) = vals;
    end
    
end
mat = temp;