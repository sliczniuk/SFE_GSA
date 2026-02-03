function cache = init_Extraction_Cache(time_step)
% INIT_EXTRACTION_CACHE Initialize cached data for extraction simulation
%   cache = init_Extraction_Cache(time_step)
%
%   This function loads parameters and builds the CasADi integrator ONCE,
%   avoiding repeated file I/O and symbolic compilation.
%
%   Input:
%       time_step - time step in minutes
%
%   Output:
%       cache - struct containing all precomputed data

addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

%% Load parameters ONCE
Parameters_table = readtable('Parameters.csv');
Parameters = num2cell(Parameters_table{:,3});

%% Geometry setup
nstages = Parameters{1};

before = 0.04;
bed = 0.92;

nstagesbefore = 1:floor(before * nstages);
nstagesbed = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed * nstages);
nstagesafter = nstagesbed(end)+1 : nstages;

bed_mask = nan(nstages, 1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed) = 1;
bed_mask(nstagesafter) = 0;

%% Dimensions
Nx = 3 * nstages + 2;
Nu = 3 + numel(Parameters);

%% Extractor geometry
r = Parameters{3};
epsi = Parameters{4};
L = Parameters{6};

L_nstages = linspace(0, L, nstages);
V_slice = (L / nstages) * pi * r^2;

V_before = V_slice * numel(nstagesbefore);
V_after = V_slice * numel(nstagesafter);
V_bed = V_slice * numel(nstagesbed);

V_before_fluid = repmat(V_before * 1 / numel(nstagesbefore), numel(nstagesbefore), 1);
V_bed_fluid = repmat(V_bed * (1 - epsi) / numel(nstagesbed), numel(nstagesbed), 1);
V_after_fluid = repmat(V_after * 1 / numel(nstagesafter), numel(nstagesafter), 1);
V_fluid = [V_before_fluid; V_bed_fluid; V_after_fluid];

L_bed_after_nstages = L_nstages(nstagesbed(1):end);
L_bed_after_nstages = L_bed_after_nstages - L_bed_after_nstages(1);
L_end = L_bed_after_nstages(end);

%% Build integrator ONCE
timeStep_in_sec = time_step * 60;
f_FP = @(x, u) modelSFE_Corr(x, u, bed_mask, timeStep_in_sec);
F_FP = buildIntegrator(f_FP, [Nx, Nu], timeStep_in_sec);

%% Store everything in cache struct
cache.Parameters = Parameters;
cache.nstages = nstages;
cache.bed_mask = bed_mask;
cache.nstagesbefore = nstagesbefore;
cache.nstagesbed = nstagesbed;
cache.nstagesafter = nstagesafter;
cache.Nx = Nx;
cache.Nu = Nu;
cache.V_bed = V_bed;
cache.epsi = epsi;
cache.V_fluid = V_fluid;
cache.L_bed_after_nstages = L_bed_after_nstages;
cache.L_end = L_end;
cache.F_FP = F_FP;
cache.time_step = time_step;

fprintf('Cache initialized: Parameters loaded, integrator built (time_step = %.2f min)\n', time_step);

end
