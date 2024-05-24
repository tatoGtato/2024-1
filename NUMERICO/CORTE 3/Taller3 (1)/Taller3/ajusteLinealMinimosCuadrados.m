function [m, b] = ajusteLinealMinimosCuadrados(X, Y)
    % Número de puntos de datos
    n = length(X);

    % Calcular las sumatorias necesarias
    sum_xy = sum(X .* Y);
    sum_x = sum(X);
    sum_y = sum(Y);
    sum_x_squared = sum(X.^2);

    % Calcular la pendiente (m) y la ordenada al origen (b)
    m = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x^2);
    b = (sum_y - m * sum_x) / n;
end
