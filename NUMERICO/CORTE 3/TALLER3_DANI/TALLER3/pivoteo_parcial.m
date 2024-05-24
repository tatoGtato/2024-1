%Para cada columna k, se selecciona la fila con el valor absoluto más grande en la columna k a partir de la fila k.
%Si esta fila no es la fila k, se intercambian las filas para colocar el pivote en la posición correcta

% x: vector de soluciones.
% A: matriz de coeficientes modificada.
% b: vector de términos independientes modificado.
function [x, A, b] = pivoteo_parcial(A, b)
    [n, ~] = size(A);  % Obtener el tamaño de la matriz A
    iter = 0; 
    x = zeros(n, 1);  % Inicializar el vector solución x con ceros

    % Fase de eliminación con pivoteo parcial
    for k = 1:n - 1
        iter = iter + 1;
        [~, max_row] = max(abs(A(k:n, k)));  % Encontrar el índice de la fila con el mayor valor absoluto en la columna k
        max_row = max_row + k - 1;  % Ajustar el índice relativo al índice absoluto

        % Intercambiar filas si es necesario
        if max_row ~= k
            A([k, max_row], :) = A([max_row, k], :);  % Intercambiar filas en A
            b([k, max_row]) = b([max_row, k]);  % Intercambiar filas en b
        end

        % Realizar la eliminación de Gauss
        for i = k + 1:n
            iter = iter + 1;
            m = A(i, k) / A(k, k);  % Calcular el multiplicador
            A(i, k:n) = A(i, k:n) - m * A(k, k:n);  % Actualizar la fila i de la matriz A
            b(i) = b(i) - m * b(k);  % Actualizar el vector b
        end
    end

    % Fase de sustitución hacia atrás
    for i = n:-1:1
        x(i) = (b(i) - A(i, i+1:end) * x(i+1:end)) / A(i, i);  % Calcular el valor de x(i)
        iter = iter + 1;
    end

    % Mostrar el número total de iteraciones
    disp(['Número total de iteraciones en Pivoteo Parcial: ', num2str(iter)]);
end
