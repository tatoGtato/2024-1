clear all 
clc

Q = @(t) 9 + 5*cos(0.4*t)^2;
c = @(t) 5*exp(-0.5*t) + 2*exp(0.15*t);
M = @(t) Q(t)*c(t);

x0 = 2;
x1 = 8;
TOL = 0.0001;

%n = 1
romberg(M, x0, x1, 1, TOL)

%n = 8
romberg(M, x0, x1, 5, TOL)

%n = 16
romberg(M, x0, x1, 8, TOL);








%% FUCIONES UTILIZADAS

%RONBERG CON TOLERANCIA
function R = romberg(f, a, b, n, tolerance)
    % Paso 1: Inicializar h y el primer elemento de R
    h = b - a;
    R = zeros(n, n); % Crear una matriz n x n para almacenar los resultados
    R(1, 1) = (h / 2) * (f(a) + f(b));

    % Paso 2: Mostrar el primer valor de R
    fprintf('R(1,1) = %.8f\n', R(1, 1));

    % Paso 3: Bucle sobre i para calcular las filas de R
    for i = 2:n
        % Paso 4: Calcular R(i, 1) usando la regla del trapecio compuesta
        sum = 0;
        for k = 1:2^(i-2)
            sum = sum + f(a + (k - 0.5) * h);
        end
        R(2, 1) = 0.5 * (R(1, 1) + h * sum);

        % Paso 5: Bucle sobre j para el cálculo de extrapolación
        for j = 2:i
            R(2, j) = R(2, j - 1) + (R(2, j - 1) - R(1, j - 1)) / (4^(j - 1) - 1);
        end

        % Paso 6: Mostrar los valores de R2,j para j = 1,2,...,i
        fprintf('R(%d,:) = ', i);
        fprintf('%.8f ', R(2, 1:i));
        fprintf('\n');

        % Paso 7: Actualizar h
        h = h / 2;

        % Paso 8: Actualizar la fila 1 de R para la siguiente iteración
        R(1, :) = R(2, :);

        % Paso 9: Verificar si se alcanzó la tolerancia
        if i > 2 && abs(R(2, i) - R(1, i - 1)) < tolerance
            fprintf('Tolerancia alcanzada. Se devuelve el valor actual.\n');
            return;
        end
    end
end
