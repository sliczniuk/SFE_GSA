function outputs = applyToEachRowOptimized(func, inputMatrix)
% APPLYTOEACHROWOPTIMIZED Applies function to each row with preallocated output
%   outputs = applyToEachRowOptimized(func, inputMatrix)
%
%   Optimized version that:
%   - Uses preallocated numeric array instead of cell array
%   - Avoids cell2mat conversion overhead
%   - Skips redundant parallel pool checks
%
%   Inputs:
%       func        - function handle @(x) that returns a scalar
%       inputMatrix - M-by-N matrix where each row is an input
%
%   Output:
%       outputs     - M-by-1 vector of outputs

    numRows = size(inputMatrix, 1);
    outputs = zeros(numRows, 1);  % Preallocate numeric array

    parfor i = 1:numRows
        outputs(i) = func(inputMatrix(i, :));
    end
end
