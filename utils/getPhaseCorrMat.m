function [phaseCorrMat] = getPhaseCorrMat(A, B)

% INPUTS
%
% A, B - nCellxN matrices of phase values, both the same size
%
% OUTPUTS
% 
% phaseCorrMat

nBins = size(A, 2);

temp = nan(nBins);

for ii = 1:nBins
    for jj = 1:nBins
        ids = ~isnan(A(:, ii)) & ~isnan(B(:, jj));
        
        temp(ii, jj) = circ_corrcc(A(ids, ii), B(ids, jj));
        
    end
end

phaseCorrMat = temp; %tril(temp.', 1) + triu(temp);
