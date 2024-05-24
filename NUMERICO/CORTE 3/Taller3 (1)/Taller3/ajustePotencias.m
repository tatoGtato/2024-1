function [a, b] = ajustePotencias(X, Y)
    % Definir la función de modelo
    modelo = @(p, x) p(1) .* x .^ p(2);

    % Estimación inicial de los parámetros
    p0 = [1, 1];

    % Ajustar los parámetros utilizando mínimos cuadrados
    parametros = lsqcurvefit(modelo, p0, X, Y);

    % Extraer los valores ajustados de los parámetros
    a = parametros(1);
    b = parametros(2);
end
