function [X, Y_aprox] = aproximacionFourier(f, a, b, N, num_puntos)
    % Coeficientes de la serie de Fourier
    a0 = (1 / (b - a)) * integral(@(x) f(x), a, b);
    an = zeros(1, N);
    bn = zeros(1, N);
    
    for n = 1:N
        an(n) = (1 / (b - a)) * integral(@(x) f(x) .* cos((n * pi * x) / (b - a)), a, b);
        bn(n) = (1 / (b - a)) * integral(@(x) f(x) .* sin((n * pi * x) / (b - a)), a, b);
    end
    
    % Aproximación de Fourier
    X = linspace(a, b, num_puntos);
    Y_aprox = zeros(1, num_puntos);
    
    for i = 1:num_puntos
        Y_aprox(i) = a0 / 2;
        for n = 1:N
            Y_aprox(i) = Y_aprox(i) + an(n) * cos((n * pi * X(i)) / (b - a)) + bn(n) * sin((n * pi * X(i)) / (b - a));
        end
    end
end
