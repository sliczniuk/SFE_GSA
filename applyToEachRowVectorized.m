function outputs = applyToEachRowVectorized(func, inputList, options)
% Vectorized version - fastest when func supports it
    arguments
        func (1,1) function_handle
        inputList
        options.BatchSize (1,1) double = 1000
        options.TryVectorize (1,1) logical = true
    end
    
    numRows = size(inputList, 1);
    
    % Try full vectorization first
    if options.TryVectorize
        try
            outputs = func(inputList);
            return;
        catch
            % Fall back to batched processing
        end
    end
    
    % Batched processing
    numBatches = ceil(numRows / options.BatchSize);
    outputs = cell(numBatches, 1);
    
    parfor i = 1:numBatches
        startIdx = (i - 1) * options.BatchSize + 1;
        endIdx = min(i * options.BatchSize, numRows);
        
        batchData = inputList(startIdx:endIdx, :);
        outputs{i} = func(batchData);
    end
    
    outputs = vertcat(outputs{:});
end