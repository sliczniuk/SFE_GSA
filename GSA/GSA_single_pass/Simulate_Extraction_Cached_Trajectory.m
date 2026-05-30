function Y = Simulate_Extraction_Cached_Trajectory(par, time_points, time_step, cache)
% SIMULATE_EXTRACTION_CACHED_TRAJECTORY Run once and return yields at requested times.
%   Y = Simulate_Extraction_Cached_Trajectory(par, time_points, time_step, cache)
%   returns a 1-by-numel(time_points) row vector. The model is integrated
%   only to max(time_points), and intermediate yields are read from the
%   simulated state trajectory.

import casadi.*

time_points = time_points(:)';
if isempty(time_points)
    error('time_points must contain at least one requested output time.');
end
if any(time_points < 0)
    error('time_points must be non-negative.');
end
if time_step <= 0
    error('time_step must be positive.');
end

step_idx = round(time_points ./ time_step);
tolerance = max(1e-10, 1e-9 * max(1, max(time_points)));
if any(abs(step_idx .* time_step - time_points) > tolerance)
    error('All time_points must align with time_step %.12g. No interpolation is applied.', time_step);
end

final_time = max(time_points);
N_Time = round(final_time / time_step);

T = par(1);
P = par(2);
F = par(3);

Parameters = cache.Parameters;
nstages = cache.nstages;
bed_mask = cache.bed_mask;
V_bed = cache.V_bed;
epsi = cache.epsi;
V_fluid = cache.V_fluid;
L_bed_after_nstages = cache.L_bed_after_nstages;
L_end = cache.L_end;
nstagesbefore = cache.nstagesbefore;
F_FP = cache.F_FP;

%% Initial conditions
m_total = 3.0;
msol_max = m_total;
mSol_ratio = 1;

mSOL_s = msol_max * mSol_ratio;
mSOL_f = msol_max * (1 - mSol_ratio);

C0solid = mSOL_s * 1e-3 / (V_bed * epsi);
Parameters{2} = C0solid;

G = @(x) -(2 * mSOL_f / L_end^2) * (x - L_end);
m_fluid = G(L_bed_after_nstages) * L_bed_after_nstages(2);
m_fluid = [zeros(1, numel(nstagesbefore)) m_fluid];
C0fluid = m_fluid * 1e-3 ./ V_fluid';

%% Operating conditions
T0homog = T + 273;
Press = P;
Flow = F;

feedTemp = T0homog * ones(1, N_Time);
feedPress = Press * ones(1, N_Time);
feedFlow = Flow * ones(1, N_Time);

uu = [feedTemp', feedPress', feedFlow'];

%% Initial state
Z = Compressibility(T0homog, feedPress(1), Parameters);
rho = rhoPB_Comp(T0homog, feedPress(1), Z, Parameters);
enthalpy_rho = rho .* SpecificEnthalpy(T0homog, feedPress(1), Z, rho, Parameters);

x0 = [C0fluid';
      C0solid * bed_mask;
      enthalpy_rho * ones(nstages, 1);
      feedPress(1);
      0];

%% Simulate once and collect requested yield samples
Parameters_init_time = [uu repmat(cell2mat(Parameters), 1, N_Time)'];
X = simulateSystem(F_FP, [], x0, Parameters_init_time);

% simulateSystem includes the initial state in column 1, so time t maps to
% integration step t/time_step plus one.
Y = X(end, step_idx + 1);

end
