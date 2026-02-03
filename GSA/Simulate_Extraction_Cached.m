function Y = Simulate_Extraction_Cached(par, final_time, time_step, cache)
% SIMULATE_EXTRACTION_CACHED Optimized extraction simulation with caching
%   Y = Simulate_Extraction_Cached(par, final_time, time_step, cache)
%
%   Inputs:
%       par        - [T, P, F] operating conditions
%       final_time - simulation end time (minutes)
%       time_step  - time step (minutes)
%       cache      - struct with precomputed data (from init_Extraction_Cache)
%
%   Output:
%       Y - yield at final_time

import casadi.*

T = par(1);
P = par(2);
F = par(3);

% Use cached parameters
Parameters = cache.Parameters;
nstages = cache.nstages;
bed_mask = cache.bed_mask;
V_bed = cache.V_bed;
epsi = cache.epsi;
V_fluid = cache.V_fluid;
L_bed_after_nstages = cache.L_bed_after_nstages;
L_end = cache.L_end;
nstagesbefore = cache.nstagesbefore;
Nx = cache.Nx;

% Use cached integrator
F_FP = cache.F_FP;

%% Time setup
timeStep_in_sec = time_step * 60;
Time_in_sec = (time_step:time_step:final_time) * 60;
N_Time = length(Time_in_sec);

%% Initial conditions (these depend on parameters, compute each time)
m_total = 3.0;
msol_max = m_total;
mSol_ratio = 1;

mSOL_s = msol_max * mSol_ratio;
mSOL_f = msol_max * (1 - mSol_ratio);

C0solid = mSOL_s * 1e-3 / (V_bed * epsi);
Parameters{2} = C0solid;

G = @(x) -(2*mSOL_f / L_end^2) * (x - L_end);
m_fluid = G(L_bed_after_nstages) * (L_bed_after_nstages(2));
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
T_0 = T0homog;
Z = Compressibility(T_0, feedPress(1), Parameters);
rho = rhoPB_Comp(T_0, feedPress(1), Z, Parameters);
enthalpy_rho = rho .* SpecificEnthalpy(T_0, feedPress(1), Z, rho, Parameters);

x0 = [C0fluid';
      C0solid * bed_mask;
      enthalpy_rho * ones(nstages, 1);
      feedPress(1);
      0];

%% Simulate
Parameters_init_time = [uu repmat(cell2mat(Parameters), 1, N_Time)'];
X = simulateSystem(F_FP, [], x0, Parameters_init_time);
Y = X(end, end);

end
