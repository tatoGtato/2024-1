function y_aprox = aproximaFourier(X, Y, grado)
    % Calcula los coeficientes de la serie de Fourier
    N = length(X);
    coeficientes = zeros(1, grado + 1);
    for k = 0:grado
        coeficientes(k + 1) = (2 / N) * sum(Y .* cos(2 * pi * k * X));
    end
    
    % Reconstruye la función aproximada utilizando la serie de Fourier
    y_aprox = zeros(size(X));
    for k = 0:grado
        y_aprox = y_aprox + coeficientes(k + 1) * cos(2 * pi * k * X);
    end
end
