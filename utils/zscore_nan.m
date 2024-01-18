function [out] = zscore_nan(inputMatrix)

% inputs
% 
% inputMatrix - Matrix to be z-scored
% 
% output
% 
% out - zscored matrix
% 
% Notes:
% Matrix is z-scored along column

out = nan(size(inputMatrix));

for ii = 1:size(inputMatrix, 2)
    out(:,ii) = (inputMatrix(:,ii)-nanmean(inputMatrix(:,ii)))./nanstd(inputMatrix(:,ii));
end




