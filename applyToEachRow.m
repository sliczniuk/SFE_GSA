function outputs = applyToEachRow(func, inputList)
% APPLYTOEACHROW Applies a function to each row of the input list
%   outputs = applyToEachRow(func, inputList)
%   - func: a function handle @(x) that accepts a row vector
%   - inputList: an M-by-N matrix or cell array of inputs
%   - outputs: cell array or numeric array of outputs

    numRows = size(inputList, 1);
    outputs = cell(numRows, 1);  % Preallocate outputs

    % Start parallel pool if not already running
    if isempty(gcp('nocreate'))
        parpool;
    end

    parfor i = 1:numRows
        if iscell(inputList)
            rowInput = inputList{i};
        else
            rowInput = inputList(i, :);
        end
        outputs{i} = func(rowInput);
    end

    % If all outputs are numeric scalars, convert to vector
    if all(cellfun(@(x) isscalar(x) && isnumeric(x), outputs))
        outputs = cell2mat(outputs);
    end
end
