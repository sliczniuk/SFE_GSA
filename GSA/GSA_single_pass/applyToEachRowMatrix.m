function outputs = applyToEachRowMatrix(func, inputMatrix)
% APPLYTOEACHROWMATRIX Apply a vector-output function to each input row.
%   outputs = applyToEachRowMatrix(func, inputMatrix) returns an
%   numRows-by-numOutputs matrix.

numRows = size(inputMatrix, 1);
if numRows == 0
    outputs = [];
    return;
end

firstOutput = func(inputMatrix(1, :));
firstOutput = firstOutput(:)';

outputs = zeros(numRows, numel(firstOutput));
outputs(1, :) = firstOutput;

if numRows > 1
    parfor i = 2:numRows
        rowOutput = func(inputMatrix(i, :));
        outputs(i, :) = rowOutput(:)';
    end
end

end
