% Para cada columna k, se encuentra la fila con el mayor ratio escalado (valor absoluto del elemento dividido por el valor máximo absoluto de la fila).
% Si la fila con el mayor ratio escalado no es la fila actual k, se intercambian las filas para colocar el pivote en la posición correcta.
% Se realiza la eliminación de Gauss para hacer ceros debajo del pivote.
function [x, A, b] = pivoteo_escalado(A, b)
    % Obtiene el tamaño de la matriz A 
    [n, ~] = size(A);
    
    % Inicializa el vector de escalamiento s
    s = zeros(n, 1);
    iter = 0;
    
    % Inicializa el vector solución x con ceros
    x = zeros(n, 1);

    % Calcula el vector de escalamiento s
    for i = 1:n
        s(i) = max(abs(A(i, :)));  % Encuentra el valor máximo absoluto de cada fila
        iter = iter + 1; 
    end

    % Fase de eliminación con pivoteo escalado
    for k = 1:n - 1
        iter = iter + 1;
        max_ratio = 0;  % Inicializa el máximo ratio
        pivot_row = k;  % Inicializa la fila pivote

        % Encuentra la fila con el mayor ratio escalado
        for i = k:n
            iter = iter + 1;
            ratio = abs(A(i, k)) / s(i);  % Calcula el ratio escalado
            if ratio > max_ratio
                max_ratio = ratio;
                pivot_row = i;  % Actualiza la fila pivote
            end
        end

        % Intercambia la fila actual con la fila pivote si es necesario
        if pivot_row ~= k
            A([k, pivot_row], :) = A([pivot_row, k], :);  % Intercambia filas en A
            b([k, pivot_row]) = b([pivot_row, k]);  % Intercambia filas en b
            s([k, pivot_row]) = s([pivot_row, k]);  % Intercambia filas en s
        end

        % Realiza la eliminación de Gauss
        for i = k + 1:n
            m = A(i, k) / A(k, k);  % Calcula el multiplicador
            A(i, k:n) = A(i, k:n) - m * A(k, k:n);  % Actualiza la fila i de la matriz A
            b(i) = b(i) - m * b(k);  % Actualiza el vector b
            iter = iter + 1;
        end
    end

    % Fase de sustitución hacia atrás
    for i = n:-1:1
        x(i) = (b(i) - A(i, i+1:n) * x(i+1:n)) / A(i, i);  % Calcula el valor de x(i)
        iter = iter + 1;
    end

    disp(['Número total de iteraciones en Pivoteo Escalado: ', num2str(iter)]);
end
