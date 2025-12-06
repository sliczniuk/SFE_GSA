function outputs = applyToEachRowGPU_arrayfun(time, inputList)
    % inputList: M×N numeric
    % time: parameter you were capturing before

    [numRows, ~] = size(inputList);

    % Move full input to GPU
    inputGPU = gpuArray(inputList);

    % Row indices on GPU
    idxGPU = gpuArray.colon(1, numRows)';

    % Call Simulate_Extraction directly from arrayfun kernel
    outputsGPU = arrayfun(@(i) kernel(i, inputGPU, time), idxGPU);

    outputs = gather(outputsGPU);
end