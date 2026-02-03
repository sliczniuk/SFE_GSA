function out = kernel(i, inputGPU, time)
    % Slice the i-th row and call the top-level function
    YY = inputGPU(i, :);
    out = Simulate_Extraction(YY, time, 1);
end