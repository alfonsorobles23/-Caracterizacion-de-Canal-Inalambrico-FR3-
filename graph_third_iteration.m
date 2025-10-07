% SCRIPT TO OPTIMIZE TWO PATH LOSS MODELS: ONE FOR LOS AND ONE FOR NLOS
% USING A FREE SPACE REFERENCE
clc;
clear;
close all;

% --- 1. System Configuration and Data Loading ---
fprintf('Loading measurement data for selected heights...\n');
f = 18e9; % Signal frequency in Hz (18 GHz)
c = 3e8;  % Speed of light in m/s
lambda = c / f; % Wavelength

% Definition of heights and corresponding files (only the selected ones)
alturas = [0.61, 1.30, 1.91]; % Receiver antenna heights in m
archivos = {'resultados_metodo_lee061.mat', 'resultados_metodo_lee130.mat', 'resultados_metodo_lee191.mat'};

% --- Pre-allocation of memory ---
temp_data = load(archivos{1});
n_los_points = length(temp_data.distancias_los);
n_nlos_points = length(temp_data.distancias_nlos);
total_files = length(archivos);
dist_total_los = zeros(n_los_points * total_files, 1);
pl_medido_los = zeros(n_los_points * total_files, 1);
dist_total_nlos = zeros(n_nlos_points * total_files, 1);
pl_medido_nlos = zeros(n_nlos_points * total_files, 1);

los_idx = 1;
nlos_idx = 1;

for i = 1:total_files
    try
        datos = load(archivos{i});
        fprintf('File "%s" loaded for h_r = %.2f m.\n', archivos{i}, alturas(i));

        los_end_idx = los_idx + length(datos.distancias_los) - 1;
        nlos_end_idx = nlos_idx + length(datos.distancias_nlos) - 1;

        dist_total_los(los_idx:los_end_idx) = datos.distancias_los;
        pl_medido_los(los_idx:los_end_idx) = datos.pl_lee_los;
        dist_total_nlos(nlos_idx:nlos_end_idx) = datos.distancias_nlos;
        pl_medido_nlos(nlos_idx:nlos_end_idx) = datos.pl_lee_nlos;

        los_idx = los_end_idx + 1;
        nlos_idx = nlos_end_idx + 1;
    catch
        warning('Could not find file "%s". Skipping to next height.', archivos{i});
    end
end
dist_total_los = dist_total_los(1:los_idx-1);
pl_medido_los = pl_medido_los(1:los_idx-1);
dist_total_nlos = dist_total_nlos(1:nlos_idx-1);
pl_medido_nlos = pl_medido_nlos(1:nlos_idx-1);

dist_total_global = [dist_total_los; dist_total_nlos];

% --- 2. Optimization of two separate CI models (free space reference) ---
fprintf('\nStarting optimization for separate Path Loss models (free space reference)...\n');

% Define the reference distance (d0) and calculate the free space reference PL
d0 = 3.15; % Reference distance in meters
PL_d0_fs = 20 * log10(4 * pi * d0 / lambda);

% 2.1 Optimize 'n_los' with LOS data only
fprintf('Optimizing path loss exponent (n_los) with LOS data...\n');
objective_function_los = @(n) ...
    sqrt(mean((pl_medido_los - (PL_d0_fs + 10 * n * log10(dist_total_los / d0))).^2));
n_los_opt = fminbnd(objective_function_los, 2, 6);

% 2.2 Optimize 'n_nlos' with NLOS data only
fprintf('Optimizing path loss exponent (n_nlos) with NLOS data...\n');
objective_function_nlos = @(n) ...
    sqrt(mean((pl_medido_nlos - (PL_d0_fs + 10 * n * log10(dist_total_nlos / d0))).^2));
n_nlos_opt = fminbnd(objective_function_nlos, 2, 6);

fprintf('\n--- Optimized Parameters for Separate Models ---\n');
fprintf('Free space reference PL (PL_d0_fs): %.2f dB\n', PL_d0_fs);
fprintf('LOS Path Loss Exponent (n_los): %.2f\n', n_los_opt);
fprintf('NLOS Path Loss Exponent (n_nlos): %.2f\n', n_nlos_opt);

% --- PRINT FINAL MODELS WITH VALUES ---
fprintf('\n--- Final Path Loss Models ---\n');
fprintf('CI Model for LOS: PL = %.2f + 10 * %.2f * log10(d / %.2f)\n', PL_d0_fs, n_los_opt, d0);
fprintf('CI Model for NLOS: PL = %.2f + 10 * %.2f * log10(d / %.2f)\n', PL_d0_fs, n_nlos_opt, d0);
% --------------------------------------------------

% --- 3. Generate Predictions and Calculate Final RMSE ---
pl_modelo_los = PL_d0_fs + 10 * n_los_opt * log10(dist_total_los / d0);
pl_modelo_nlos = PL_d0_fs + 10 * n_nlos_opt * log10(dist_total_nlos / d0);
rmse_los = sqrt(mean((pl_medido_los - pl_modelo_los).^2));
rmse_nlos = sqrt(mean((pl_medido_nlos - pl_modelo_nlos).^2));

pl_pred_total = [pl_modelo_los; pl_modelo_nlos];
pl_medido_total_sorted = [pl_medido_los; pl_medido_nlos];
rmse_global_final = sqrt(mean((pl_medido_total_sorted - pl_pred_total).^2));

fprintf('\n----------- FINAL RMSE TABLE [dB] -----------\n');
fprintf('Separate Model     |  LOS    |  NLOS  |  Both   \n');
fprintf('-------------------|---------|--------|---------\n');
fprintf('Optimized CI       |  %5.2f  |  %5.2f |  %5.2f\n', rmse_los, rmse_nlos, rmse_global_final);
fprintf('--------------------------------------------------\n');

% --- 4. Final Graph Generation (Linear Scale) ---
fprintf('Generating linear scale graph...\n');
figure('Name', 'Optimized CI Model vs. Measurement (Linear Scale)');
hold on;
font_format = {'FontName', 'Times New Roman', 'FontSize', 14};
poster_line_width = 3.5;
poster_marker_size = 4;
colors_medido = {[0.5 0.2 0.8], [0.8 0.5 0.2], [0.1 0.7 0.5]};

for i = 1:length(alturas)
    datos = load(archivos{i});
    dist_temp_los = datos.distancias_los;
    pl_temp_los = datos.pl_lee_los;
    plot(dist_temp_los, pl_temp_los, 'o', ...
         'MarkerSize', poster_marker_size, ...
         'MarkerEdgeColor', colors_medido{i}, ...
         'MarkerFaceColor', colors_medido{i}, ...
         'DisplayName', sprintf('LOS Data h_r=%.2fm', alturas(i)));
    
    dist_temp_nlos = datos.distancias_nlos;
    pl_temp_nlos = datos.pl_lee_nlos;
    plot(dist_temp_nlos, pl_temp_nlos, 's', ...
         'MarkerSize', poster_marker_size, ...
         'MarkerEdgeColor', colors_medido{i}, ...
         'MarkerFaceColor', colors_medido{i}, ...
         'DisplayName', sprintf('NLOS Data h_r=%.2fm', alturas(i)));
end

los_range = [min(dist_total_los), max(dist_total_los)];
nlos_range = [min(dist_total_nlos), max(dist_total_nlos)];

manhattan_dist_los_model = linspace(los_range(1), los_range(2), 200);
pl_modelo_global_los = PL_d0_fs + 10 * n_los_opt * log10(manhattan_dist_los_model / d0);
plot(manhattan_dist_los_model, pl_modelo_global_los, 'k-', 'LineWidth', poster_line_width, 'DisplayName', 'Optimized CI Model (LOS)');

manhattan_dist_nlos_model = linspace(nlos_range(1), nlos_range(2), 200);
pl_modelo_global_nlos = PL_d0_fs + 10 * n_nlos_opt * log10(manhattan_dist_nlos_model / d0);
plot(manhattan_dist_nlos_model, pl_modelo_global_nlos, 'r-', 'LineWidth', poster_line_width, 'DisplayName', 'Optimized CI Model (NLOS)');

dist_friis = linspace(min(dist_total_global), max(dist_total_global), 200);
pl_friis = 20 * log10(4 * pi * dist_friis / lambda);
plot(dist_friis, pl_friis, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', poster_line_width - 1, 'DisplayName', 'Free Space (Friis)');

hold off;
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.4;
ax.Layer = 'top';

% Modified legend position and font size
leg = legend('Location', 'northwest', 'NumColumns', 2);
leg.FontSize = 12; % Increased font size
leg.Position = [0.08, 0.76, 0.35, 0.18]; % Adjusted position for a larger box

% Set titles and labels in English
title('Path Loss vs. Distance: Global Model vs. Measurements', 'FontSize', 16);
xlabel('Manhattan Distance [m]');
ylabel('Path Loss [dB]');

xlim([min(dist_total_global), max(dist_total_global)]);
ylim([40 160]);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);

fprintf('Process finished.\n');