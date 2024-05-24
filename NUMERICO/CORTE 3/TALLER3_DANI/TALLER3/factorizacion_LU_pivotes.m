% pivotes parciales (también conocidos como pivotes de intercambio) para hacer la factorización LU más robusta y capaz
% de manejar matrices que no son estrictamente diagonalmente dominantes o que tienen ceros en la diagonal. 
% La inclusión de pivotes parciales asegura que en cada paso del proceso, se elija el pivote más grande en magnitud 
% de los disponibles en la columna actual, lo que ayuda a evitar la división por cero y reduce los errores de redondeo.
function [L, U, P] = factorizacion_LU_pivotes(A)
    [n, ~] = size(A);
    L = eye(n);
    U = zeros(n);
    P = eye(n); % Matriz de permutación inicializada como identidad

    for k = 1:n
        % Buscar el pivote más grande en magnitud en la columna k
        [~, p] = max(abs(A(k:end, k)));
        p = p + k - 1; % Índice global del pivote

        % Intercambiar filas en A, L y P si es necesario
        if p ~= k
            A([k, p], :) = A([p, k], :);
            L([k, p], :) = L([p, k], :);
            P([k, p], :) = P([p, k], :);
        end

        % Continuar con la factorización LU
        for j = k:n
            U(k, j) = A(k, j) - L(k, 1:k-1) * U(1:k-1, j);
        end

        for i = k+1:n
            L(i, k) = (A(i, k) - L(i, 1:k-1) * U(1:k-1, k)) / U(k, k);
        end
    end
    
    % Aplicar la matriz de permutación a la matriz A original
    A_perm = P * A;
    
    % Verificar que la factorización sea correcta: A = L * U
    % L * U debería ser igual a A_perm (A con filas permutadas)
    LU = L * U;
    disp('Verificación de la factorización LU con pivotes:');
    disp('A_perm - LU debería ser una matriz de ceros:');
    disp(A_perm - LU);
end