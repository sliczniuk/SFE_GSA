% Sobol Sensitivity Analysis in MATLAB
% General-purpose script applicable to any deterministic model

close all; clear; clc;

%% USER INPUTS

alpha      = 0.15;
case_name  = 2;
FONT       = 10;

% Number of input parameters (dimensions)
d = 3; % Now extended to 3 parameters

% Number of Monte Carlo samples (recommended: thousands for accuracy)
N = 1000;

% Define min and max values for parameters (uniform distributions)
param_min = [30+273, 100, 3.33e-5]; % Min values for parameters 1, 2, and 3
param_max = [40+273, 200, 6.66e-5]; % Max values for parameters 1, 2, and 3

NAMES     = {'Temperature', 'Pressure', 'Mass flowrate'};

% Model function handle (USER: modify this to your own model)

switch case_name
    case 1
        model_func = @(X) Compute_Gamma_for_GSA(X);
        case_name  = '$\Upsilon$';
        save_name  = 'Upsilon';
    case 2
        model_func = @(X) Compute_Diffusion_for_GSA(X);
        case_name  = '$D_i^R$';
        save_name  = 'DiR';
    case 3
        model_func = @(X) Compute_kinetic_for_GSA(X, 1);
        case_name  = '$D_i$';
        save_name  = 'Di';
    otherwise
        warning('Wrong case')
end

%% STEP 1: Generate Sobol samples

A = rand(N, d);
B = rand(N, d);
for i = 1:d
    A(:,i) = param_min(i) + (param_max(i) - param_min(i)) .* A(:,i);
    B(:,i) = param_min(i) + (param_max(i) - param_min(i)) .* B(:,i);
end

YA = model_func(A);
YB = model_func(B);

Y = [YA; YB];
VY = var(Y);

%% STEP 2: Compute First-order and Total-order Sobol indices, store YABi

first_order = zeros(d,1);
total_order = zeros(d,1);
YAB = zeros(N,d);

for i = 1:d
    ABi = A;
    ABi(:, i) = B(:, i);
    YAB(:,i) = model_func(ABi);
    first_order(i) = (1/N)*sum(YB .* (YAB(:,i) - YA)) / VY;
    total_order(i) = (1/(2*N))*sum((YA - YAB(:,i)).^2) / VY;
end

%% STEP 2b: Compute Second-order Sobol indices

second_order = zeros(d,d);

for i = 1:d-1
    for j = i+1:d
        ABij = A;
        ABij(:, i) = B(:, i);
        ABij(:, j) = B(:, j);
        YABij = model_func(ABij);
        second_order(i,j) = (1/N)*sum(YB .* (YABij - YAB(:,i) - YAB(:,j) + YA)) / VY;
    end
end

%% STEP 3: Display Results

fprintf('\nSobol Sensitivity Indices:\n');
for i = 1:d
    fprintf('%s: First-order = %.4f, Total-order = %.4f\n', ...
        NAMES{i}, first_order(i), total_order(i));
end

fprintf('\nSecond-order Sobol Indices (upper triangle, i<j):\n');
for i = 1:d-1
    for j = i+1:d
        fprintf('%s & %s: Second-order = %.4f\n', NAMES{i}, NAMES{j}, second_order(i,j));
    end
end

%% STEP 4: Bar plot for First, Second, and Total-order indices
%{
figure;
bar([first_order, total_order]);
hold on;
for i = 1:d-1
    for j = i+1:d
        plot([i, j], [second_order(i,j), second_order(i,j)], 'ko-', 'LineWidth',2);
    end
end
%xlabel('Parameter');
xticklabels(NAMES)
ylabel('Sobol Sensitivity Index');
%title('Sobol Sensitivity Indices (1st, 2nd, Total)');
lgd = legend('First-order', 'Total-order', 'Second-order pairs', Location='northwest');
lgd.FontSize = FONT; lgd.Box='off';

grid on; hold off;
set(gca,'FontSize',FONT)
exportgraphics(figure(1), ['GSA_BAR_U',save_name,'.png'], "Resolution",300); close all;

%% STEP 5: Scatter plots for parameter pairs vs. Output

figure;
%subplot(1,3,1);
h=scatter3(A(:,1), A(:,2), YA, 10, YA, 'filled');
set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
set(gca,'FontSize',FONT)
xlabel('Tempertaure [K]'); ylabel('Pressure [bar]'); zlabel(case_name);
%title('Parameters 1 & 2 vs. Output'); 
grid on; colorbar;
%view(-18, 20)
view(145, 25)
axis square
exportgraphics(figure(1), ['GSA_Scatter_1_',save_name,'.png'], "Resolution",300); close all;

%subplot(1,3,2);
figure;
h=scatter3(A(:,1), A(:,3), YA, 10, YA, 'filled');
set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
set(gca,'FontSize',FONT)
xlabel('Temperature [K]'); ylabel('Mass flowrate [kg/s]'); zlabel(case_name);
%view(-18, 20)
view(145, 25)
axis square
%title('Parameters 1 & 3 vs. Output'); 
grid on; colorbar;
exportgraphics(figure(1), ['GSA_Scatter_2_',save_name,'.png'], "Resolution",300); close all;

%subplot(1,3,3);
figure;
h = scatter3(A(:,2), A(:,3), YA, 10, YA, 'filled');
set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
set(gca,'FontSize',FONT)
xlabel('Pressure [bar]'); ylabel('Mass flowrate [kg/s]'); zlabel(case_name);
%title('Parameters 2 & 3 vs. Output'); 
grid on; colorbar;
%view(-18, 20)
view(145, 25)
axis square
exportgraphics(figure(1), ['GSA_Scatter_3_',save_name,'.png'], "Resolution",300); close all;

%% STEP 6: Inputs vs Outputs

for i = 1:d
    %subplot(1,d,i);
    scatter(A(:,i), YA, 10, YA, 'filled', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha); colorbar;
    xlabel(sprintf('%s', NAMES{i}));
    ylabel(case_name);
    %title(sprintf('Parameter %d vs. Output', i));
    grid off; axis square;
    set(gca,'FontSize',FONT)
    exportgraphics(figure(1), ['GSA_Scatter_',num2str(3+i),'_',save_name,'.png'], "Resolution",300); close all;
end

%% STEP 7: Compute and Plot Output Distribution

% Compute output samples (already computed YA and YB can be combined)
all_outputs = [YA; YB];

% Plot histogram of output distribution
figure;
histogram(all_outputs, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.6, 0.8]);
xlabel(case_name);
ylabel('Probability Density');
%title('Output Probability Distribution');
grid off;

% Compute statistics (mean and standard deviation)
output_mean = mean(all_outputs);
output_std = std(all_outputs);

% Display statistics
fprintf('\nOutput Statistics:\n');
fprintf('Mean of output: %.4f\n', output_mean);
fprintf('Standard deviation of output: %.4f\n', output_std);
set(gca,'FontSize',FONT)
exportgraphics(figure(1), ['GSA_Distribution_',save_name,'.png'], "Resolution",300); close all;
%}