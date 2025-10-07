% SCRIPT PARA OPTIMIZAR EL MODELO CI UTILIZANDO EL MÉTODO DEL PAPER
clc;
clear;
close all;
% --- 1. Configuración del Sistema y Carga de Datos ---
fprintf('Cargando datos de mediciones de las alturas seleccionadas...\n');
f = 18e9; % Frecuencia de la señal en Hz (18 GHz)
c = 3e8;  % Velocidad de la luz en m/s
lambda = c / f; % Longitud de onda
% Definición de las alturas y los archivos correspondientes (solo las seleccionadas)
alturas = [0.61, 1.30, 1.91]; % Alturas de la antena receptora en m
archivos = {'resultados_metodo_lee061.mat', 'resultados_metodo_lee130.mat', 'resultados_metodo_lee191.mat'};
% --- SOLUCIÓN PARA LOS WARNINGS: Pre-asignación de memoria ---
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
        fprintf('Archivo "%s" cargado para h_r = %.2f m.\n', archivos{i}, alturas(i));
        los_end_idx = los_idx + length(datos.distancias_los) - 1;
        nlos_end_idx = nlos_idx + length(datos.distancias_nlos) - 1;
        dist_total_los(los_idx:los_end_idx) = datos.distancias_los;
        pl_medido_los(los_idx:los_end_idx) = datos.pl_lee_los;
        dist_total_nlos(nlos_idx:nlos_end_idx) = datos.distancias_nlos;
        pl_medido_nlos(nlos_idx:nlos_end_idx) = datos.pl_lee_nlos;
        los_idx = los_end_idx + 1;
        nlos_idx = nlos_end_idx + 1;
    catch
        warning('No se pudo encontrar el archivo "%s". Saltando a la siguiente altura.', archivos{i});
    end
end
dist_total_los = dist_total_los(1:los_idx-1);
pl_medido_los = pl_medido_los(1:los_idx-1);
dist_total_nlos = dist_total_nlos(1:nlos_idx-1);
pl_medido_nlos = pl_medido_nlos(1:nlos_idx-1);
dist_total_global = [dist_total_los; dist_total_nlos];
% --- 2. Ajuste de modelos por mínimos cuadrados según el paper ---
fprintf('\nIniciando ajuste de modelos segun la metodologia del paper...\n');
d0 = 3.15; % Distancia de referencia
% --- 2.1 Ajuste para los datos de LOS ---
x_los = 10 * log10(dist_total_los / d0);
y_los = pl_medido_los;
% Ajuste lineal: y = m*x + b
coefficients_los = polyfit(x_los, y_los, 1);
n_los_opt = coefficients_los(1);
A_los = coefficients_los(2);
% --- 2.2 Ajuste para los datos de NLOS ---
x_nlos = 10 * log10(dist_total_nlos / d0);
y_nlos = pl_medido_nlos;
% Ajuste lineal: y = m*x + b
coefficients_nlos = polyfit(x_nlos, y_nlos, 1);
n_nlos_opt = coefficients_nlos(1);
A_nlos = coefficients_nlos(2);
fprintf('\n--- Parámetros Optimizados para modelos separados ---\n');
fprintf('Exponente de Pérdida LOS (n_los): %.2f\n', n_los_opt);
fprintf('Exponente de Pérdida NLOS (n_nlos): %.2f\n', n_nlos_opt);
fprintf('Intercepto LOS (A_los): %.2f dB\n', A_los);
fprintf('Intercepto NLOS (A_nlos): %.2f dB\n', A_nlos);
% --- IMPRIMIR LOS MODELOS FINALES CON VALORES ---
fprintf('\n--- Modelos Finales de Path Loss ---\n');
fprintf('Modelo CI para LOS: PL = %.2f + 10 * %.2f * log10(d / %.2f)\n', A_los, n_los_opt, d0);
fprintf('Modelo CI para NLOS: PL = %.2f + 10 * %.2f * log10(d / %.2f)\n', A_nlos, n_nlos_opt, d0);
% --------------------------------------------------
% --- 3. Generar Predicciones y Calcular RMSE Final ---
pl_modelo_los = A_los + 10 * n_los_opt * log10(dist_total_los / d0);
pl_modelo_nlos = A_nlos + 10 * n_nlos_opt * log10(dist_total_nlos / d0);
rmse_los = sqrt(mean((pl_medido_los - pl_modelo_los).^2));
rmse_nlos = sqrt(mean((pl_medido_nlos - pl_modelo_nlos).^2));
pl_pred_total = [pl_modelo_los; pl_modelo_nlos];
pl_medido_total = [pl_medido_los; pl_medido_nlos];
rmse_global_final = sqrt(mean((pl_medido_total - pl_pred_total).^2));
fprintf('\n----------- TABLA DE RMSE FINAL [dB] -----------\n');
fprintf('Modelo Separado    |  LOS    |  NLOS  |  Ambos  \n');
fprintf('-------------------|---------|--------|---------\n');
fprintf('CI Optimizado      |  %5.2f  |  %5.2f |  %5.2f\n', rmse_los, rmse_nlos, rmse_global_final);
fprintf('--------------------------------------------------\n');
% --- 4. Generación del Gráfico Final (Escala Lineal) ---
fprintf('Generando gráfico en escala lineal...\n');
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
         'DisplayName', sprintf('Datos LOS h_r=%.2fm', alturas(i)));
    dist_temp_nlos = datos.distancias_nlos;
    pl_temp_nlos = datos.pl_lee_nlos;
    plot(dist_temp_nlos, pl_temp_nlos, 's', ...
         'MarkerSize', poster_marker_size, ...
         'MarkerEdgeColor', colors_medido{i}, ...
         'MarkerFaceColor', colors_medido{i}, ...
         'DisplayName', sprintf('Datos NLOS h_r=%.2fm', alturas(i)));
end
los_range = [min(dist_total_los), max(dist_total_los)];
nlos_range = [min(dist_total_nlos), max(dist_total_nlos)];
manhattan_dist_los_model = linspace(los_range(1), los_range(2), 200);
pl_modelo_global_los = A_los + 10 * n_los_opt * log10(manhattan_dist_los_model / d0);
plot(manhattan_dist_los_model, pl_modelo_global_los, 'k-', 'LineWidth', poster_line_width, 'DisplayName', 'Modelo CI Optimizado (LOS)');
manhattan_dist_nlos_model = linspace(nlos_range(1), nlos_range(2), 200);
pl_modelo_global_nlos = A_nlos + 10 * n_nlos_opt * log10(manhattan_dist_nlos_model / d0);
plot(manhattan_dist_nlos_model, pl_modelo_global_nlos, 'r-', 'LineWidth', poster_line_width, 'DisplayName', 'Modelo CI Optimizado (NLOS)');
dist_friis = linspace(min(dist_total_global), max(dist_total_global), 200);
pl_friis = 20 * log10(4 * pi * dist_friis / lambda);
plot(dist_friis, pl_friis, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', poster_line_width - 1, 'DisplayName', 'Espacio Libre (Friis)');
hold off;
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.4;
ax.Layer = 'top';
leg = legend('Location', 'northwest', 'NumColumns', 2);
leg.Position = [0.08, 0.78, 0.3, 0.15];
title('Path Loss vs. Distancia: Modelo Global vs. Mediciones', 'FontSize', 16);
xlabel('Distancia Manhattan [m]');
ylabel('Path Loss [dB]');
xlim([min(dist_total_global), max(dist_total_global)]);
ylim([40 160]);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
set(leg, 'FontSize', 10);
fprintf('Proceso finalizado.\n');