% SCRIPT TO OPTIMIZE THE CI MODEL USING THE METHOD FROM THE PAPER
clc;
clear;
close all;

% --- 1. System Configuration and Data Loading ---
fprintf('Loading measurement data for the selected heights...\n');
f = 18e9; % Signal frequency in Hz (18 GHz)
c = 3e8;  % Speed of light in m/s
lambda = c / f; % Wavelength

% Heights of the receiver antenna and corresponding files
alturas = [0.61, 1.30, 1.91]; % Receiver antenna heights in m
archivos = {'resultados_metodo_lee061.mat', ...
            'resultados_metodo_lee130.mat', ...
            'resultados_metodo_lee191.mat'};

% Pre-allocation of memory
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
        warning('File "%s" not found. Skipping to next height.', archivos{i});
    end
end
dist_total_los = dist_total_los(1:los_idx-1);
pl_medido_los = pl_medido_los(1:los_idx-1);
dist_total_nlos = dist_total_nlos(1:nlos_idx-1);
pl_medido_nlos = pl_medido_nlos(1:nlos_idx-1);
dist_total_global = [dist_total_los; dist_total_nlos];

% --- 2. CI Model Fitting According to the Paper ---
fprintf('\nStarting CI model fitting with reference at d0 = 3.15 m...\n');
d0 = 3.15; % Reference distance in meters

% FSPL at reference distance
fspl_d0 = 20*log10(4*pi*d0*f/c);

% --- 2.1 LOS Adjustment ---
A_los = pl_medido_los - fspl_d0;
D_los = 10*log10(dist_total_los / d0);
n_los_opt = sum(D_los .* A_los) / sum(D_los.^2);

% Deterministic CI model (no Xsigma)
pl_modelo_los = fspl_d0 + 10*n_los_opt*log10(dist_total_los / d0);

% Shadow fading (residuals) and std
Xsigma_los = pl_medido_los - pl_modelo_los;
sigma_los = sqrt(mean(Xsigma_los.^2));

% --- 2.2 NLOS Adjustment ---
A_nlos = pl_medido_nlos - fspl_d0;
D_nlos = 10*log10(dist_total_nlos / d0);
n_nlos_opt = sum(D_nlos .* A_nlos) / sum(D_nlos.^2);

% Deterministic CI model (no Xsigma)
pl_modelo_nlos = fspl_d0 + 10*n_nlos_opt*log10(dist_total_nlos / d0);

% Shadow fading (residuals) and std
Xsigma_nlos = pl_medido_nlos - pl_modelo_nlos;
sigma_nlos = sqrt(mean(Xsigma_nlos.^2));

% --- Report Parameters ---
fprintf('\n--- Optimized CI Parameters ---\n');
fprintf('Path Loss Exponent LOS (n_los): %.2f\n', n_los_opt);
fprintf('Path Loss Exponent NLOS (n_nlos): %.2f\n', n_nlos_opt);
fprintf('FSPL(f, d0=%.2f m): %.2f dB\n', d0, fspl_d0);
fprintf('\n--- Shadow Fading Parameters ---\n');
fprintf('LOS shadow fading std (sigma_los): %.2f dB\n', sigma_los);
fprintf('NLOS shadow fading std (sigma_nlos): %.2f dB\n', sigma_nlos);
fprintf('\n--- Showing some Xsigma values ---\n');
disp('First 5 LOS residuals (Xsigma_los):');
disp(Xsigma_los(1:min(5,end)));
disp('First 5 NLOS residuals (Xsigma_nlos):');
disp(Xsigma_nlos(1:min(5,end)));
fprintf('\n--- Final Path Loss Models (Deterministic CI + Residuals separately) ---\n');
fprintf('CI LOS:  PL = %.2f + 10 * %.2f * log10(d/%.2f) + Xsigma\n', fspl_d0, n_los_opt, d0);
fprintf('CI NLOS: PL = %.2f + 10 * %.2f * log10(d/%.2f) + Xsigma\n', fspl_d0, n_nlos_opt, d0);

% --- 3. RMSE Calculation (deterministic model vs measurements) ---
rmse_los = sqrt(mean((pl_medido_los - pl_modelo_los).^2));
rmse_nlos = sqrt(mean((pl_medido_nlos - pl_modelo_nlos).^2));
pl_pred_total = [pl_modelo_los; pl_modelo_nlos];
pl_medido_total = [pl_medido_los; pl_medido_nlos];
rmse_global_final = sqrt(mean((pl_medido_total - pl_pred_total).^2));
fprintf('\n----------- FINAL RMSE TABLE [dB] -----------\n');
fprintf('CI Deterministic  |  LOS    |  NLOS  |  Both  \n');
fprintf('------------------|---------|--------|---------\n');
fprintf('CI Optimized      |  %5.2f  |  %5.2f |  %5.2f\n', rmse_los, rmse_nlos, rmse_global_final);
fprintf('---------------------------------------------\n');

% --- 4. Final Plot (English, Times New Roman) ---
fprintf('Generating plot in linear scale...\n');
figure('Name', 'CI Model vs. Measurements (Linear Scale)');
hold on;
poster_line_width = 3.5;
poster_marker_size = 3;
colors_medido = {[0.5 0.2 0.8], [0.8 0.5 0.2], [0.1 0.7 0.5]};

for i = 1:length(alturas)
    datos = load(archivos{i});
    dist_temp_los = datos.distancias_los;
    pl_temp_los = datos.pl_lee_los;
    
    % Plot LOS data with a filled circle marker
    plot(dist_temp_los, pl_temp_los, 'o', ...
         'MarkerSize', poster_marker_size, ...
         'MarkerEdgeColor', colors_medido{i}, ...
         'MarkerFaceColor', colors_medido{i}, ...
         'DisplayName', sprintf('LOS Data h_r=%.2f m', alturas(i)));

    dist_temp_nlos = datos.distancias_nlos;
    pl_temp_nlos = datos.pl_lee_nlos;

    % Plot NLOS data with a transparent circle marker
    plot(dist_temp_nlos, pl_temp_nlos, 'o', ...
         'MarkerSize', poster_marker_size, ...
         'MarkerEdgeColor', colors_medido{i}, ...
         'MarkerFaceColor', 'none', ...
         'DisplayName', sprintf('NLOS Data h_r=%.2fm', alturas(i)));
end

% LOS CI Curve (deterministic)
los_range = [min(dist_total_los), max(dist_total_los)];
manhattan_dist_los_model = linspace(los_range(1), los_range(2), 200);
pl_modelo_global_los = fspl_d0 + 10*n_los_opt*log10(manhattan_dist_los_model / d0);
plot(manhattan_dist_los_model, pl_modelo_global_los, 'k-', ...
     'LineWidth', poster_line_width, 'DisplayName', 'CI Model: LOS');

% NLOS CI Curve (deterministic)
nlos_range = [min(dist_total_nlos), max(dist_total_nlos)];
manhattan_dist_nlos_model = linspace(nlos_range(1), nlos_range(2), 200);
pl_modelo_global_nlos = fspl_d0 + 10*n_nlos_opt*log10(manhattan_dist_nlos_model / d0);
plot(manhattan_dist_nlos_model, pl_modelo_global_nlos, 'r-', ...
     'LineWidth', poster_line_width, 'DisplayName', 'CI Model: NLOS');

% Free-space Friis Curve
dist_friis = linspace(min(dist_total_global), max(dist_total_global), 200);
pl_friis = 20 * log10(4 * pi * dist_friis / lambda);
plot(dist_friis, pl_friis, '--', 'Color', [0.3 0.3 0.3], ...
     'LineWidth', poster_line_width - 1, 'DisplayName', 'Free Space');

hold off;
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.4;
ax.Layer = 'top';

leg = legend('Location', 'northwest', 'NumColumns', 2);
leg.Position = [0.08, 0.78, 0.3, 0.15];
title('Path Loss vs. Distance: CI Model vs. Measurements', 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Manhattan Distance [m]', 'FontName', 'Times New Roman');
ylabel('Path Loss [dB]', 'FontName', 'Times New Roman');
xlim([min(dist_total_global), max(dist_total_global)]);
ylim([40 160]);

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
set(leg, 'FontSize', 10, 'FontName', 'Times New Roman');

fprintf('Process finished.\n');