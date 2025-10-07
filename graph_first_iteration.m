% SCRIPT PARA COMPARACIÓN DE MODELOS DE PATH LOSS CON MÚLTIPLES ALTURAS
% Este script carga datos de mediciones de path loss para cuatro alturas de
% antena, aplica el modelo de referencia CI (Close-In) y el modelo de Friis,
% y genera un único gráfico comparativo con el estilo especificado.
clc;
clear;
close all;

% --- 1. Configuración del Sistema y Carga de Datos ---
fprintf('Cargando datos de mediciones...\n');
f = 18e9; % Frecuencia de la señal en Hz (18 GHz)
c = 3e8;  % Velocidad de la luz en m/s
lambda = c / f; % Longitud de onda

% Definición de las alturas y los archivos correspondientes
alturas = [0.01, 0.61, 1.30, 1.91]; % Alturas de la antena receptora en m
archivos = {'resultados_metodo_lee001.mat', 'resultados_metodo_lee061.mat', 'resultados_metodo_lee130.mat', 'resultados_metodo_lee191.mat'};
% Definición de colores para cada altura en el gráfico
colors_medido = {[0.6 0.8 0.4], [0.5 0.2 0.8], [0.8 0.5 0.2], [0.1 0.7 0.5]}; % Nuevo color para 0.61m
colors_ci = {[0 0.4470 0.7410], [0.2 0.7 0.9], [0.4940 0.1840 0.5560], [0.8500 0.3250 0.0980]}; % Nuevo color para 0.61m

% --- 2. Parámetros del Modelo CI y Friis ---
d0 = 1; % Distancia de referencia para el modelo CI (1 metro)
% Definir el exponente de path loss para cada caso. Aquí se usa una suposición
% o un valor típico para cada ambiente, pero podrías optimizarlo.
n_ci = [2.2, 2.4, 2.5, 2.8]; % Exponentes de path loss para las 4 alturas

% --- 3. Generación del Gráfico y Cálculo de RMSE ---
figure('Name', 'Comparación de Path Loss vs Distancia para Múltiples Alturas');
hold on;
font_format = {'FontName', 'Times New Roman', 'FontSize', 20};
poster_line_width = 3.5;
poster_marker_size = 4;

rmse_table = zeros(length(alturas), 2); % Tabla para almacenar los RMSE [Altura | RMSE]
fprintf('\n--- Cálculo y Visualización ---\n');

% Bucle para procesar cada archivo (cada altura)
for i = 1:length(alturas)
    try
        datos = load(archivos{i});
        fprintf('Archivo "%s" cargado con éxito para h_r = %.2f m.\n', archivos{i}, alturas(i));

        % Extracción de datos
        dist_total = [datos.distancias_los; datos.distancias_nlos];
        pl_medido_total = [datos.pl_lee_los; datos.pl_lee_nlos];
        
        % Ordenar los datos por distancia para el gráfico
        [dist_total, sort_idx] = sort(dist_total);
        pl_medido_total = pl_medido_total(sort_idx);

        % Aplicación del modelo CI (Close-In)
        % PL_CI(d) = PL_d0 + 10 * n * log10(d / d0)
        PL_d0 = 20 * log10(4 * pi * d0 / lambda); % PL en el punto de referencia (d0)
        pl_ci_predicho = PL_d0 + 10 * n_ci(i) * log10(dist_total / d0);

        % Calcular el RMSE para el modelo CI
        rmse_ci = sqrt(mean((pl_medido_total - pl_ci_predicho).^2));
        rmse_table(i, 1) = alturas(i);
        rmse_table(i, 2) = rmse_ci;

        % Graficar los datos medidos y el modelo CI para esta altura
        plot(dist_total, pl_medido_total, 'o', ...
             'MarkerSize', poster_marker_size, ...
             'MarkerEdgeColor', colors_medido{i}, ...
             'MarkerFaceColor', colors_medido{i}, ...
             'DisplayName', sprintf('Datos Medidos h_r = %.2f m', alturas(i)));

        plot(dist_total, pl_ci_predicho, '-', ...
             'Color', colors_ci{i}, ...
             'LineWidth', poster_line_width, ...
             'DisplayName', sprintf('Modelo CI (n=%.2f) h_r = %.2f m', n_ci(i), alturas(i)));

    catch
        warning('No se pudo encontrar el archivo "%s". Saltando a la siguiente altura.', archivos{i});
    end
end

% --- 4. Graficar el Modelo de Friis ---
dist_friis = linspace(0.1, 60, 200);
pl_friis = 20 * log10(4 * pi * dist_friis / lambda);
plot(dist_friis, pl_friis, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', poster_line_width - 1, 'DisplayName', 'Espacio Libre (Friis)');

% --- 5. Ajustes Finales del Gráfico ---
hold off;
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.4;
ax.Layer = 'top';

% Ajustar la posición y el tamaño de la leyenda
leg = legend('Location', 'northwest', 'NumColumns', 2);
% Las coordenadas son [izquierda, abajo, ancho, alto] en unidades normalizadas
% Se ha ajustado para un recuadro más pequeño
leg.Position = [0.08, 0.78, 0.3, 0.15]; 

title('Path Loss vs. Distancia: Medición vs. Modelos', 'FontSize', 22); % Fuente un poco más pequeña
xlabel('Distancia [m]');
ylabel('Path Loss [dB]');
ylim([40 160]); % Ajustar el rango del eje Y para una mejor visualización

% Aplicar un tamaño de fuente más pequeño a todos los elementos
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14); % Tamaño de fuente más pequeño en general
set(leg, 'FontSize', 12); % Tamaño de fuente de la leyenda aún más pequeño

% --- 6. Presentación de Resultados y RMSE ---
fprintf('\n----------- TABLA DE RMSE POR ALTURA [dB] -----------\n');
fprintf('Altura h_r (m) | RMSE del Modelo CI\n');
fprintf('--------------------------------------------------\n');
for j = 1:size(rmse_table, 1)
    if rmse_table(j, 1) ~= 0
        fprintf('    %7.2f    |      %5.2f\n', rmse_table(j, 1), rmse_table(j, 2));
    end
end
fprintf('--------------------------------------------------\n');
fprintf('Proceso finalizado.\n');