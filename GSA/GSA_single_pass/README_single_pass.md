# GSA Single-Pass Variant

This folder contains a non-destructive refactor of the optimized GSA workflow.
The original root files, including `main_GSA_SFE_optimized.m` and
`Sen_Saltelli.m`, are not modified.

## What Changed

- `main_GSA_SFE_single_pass.m` is the new entrypoint.
- Each Saltelli sample is simulated once to `max(time_points)`.
- Intermediate yields are read from the simulated trajectory.
- Sobol indices are computed column-by-column from stored time-series outputs.
- For three parameters, the third-order index is directly estimated from the
  available Saltelli outputs rather than forced by closure.
- The model evaluation count drops from
  `(2 + d + d*(d-1)/2) * N * n_times`
  to
  `(2 + d + d*(d-1)/2) * N`.

## Files

- `main_GSA_SFE_single_pass.m` - standalone entrypoint.
- `Simulate_Extraction_Cached_Trajectory.m` - trajectory-output simulator.
- `applyToEachRowMatrix.m` - parallel row helper for vector outputs.
- `Sen_Saltelli_TimeSeries.m` - Saltelli Sobol estimator for time-series outputs.

## How to Run

From MATLAB, run:

```matlab
run(fullfile('GSA_single_pass', 'main_GSA_SFE_single_pass.m'))
```

The script adds both this folder and the parent project folder to the MATLAB
path. It temporarily changes to the parent folder while running so existing
relative dependencies such as `Parameters.csv` continue to resolve.

Results and the diary log are written inside `GSA_single_pass`.

## Assumptions

- Yield is the final state row, `X(end,:)`.
- `time_points` must align with the selected `time_step`.
- No interpolation is applied.
