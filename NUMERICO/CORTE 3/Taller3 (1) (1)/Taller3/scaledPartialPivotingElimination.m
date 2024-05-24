function x = scaledPartialPivotingElimination(A)
    [n, ~] = size(A); % Obtener el número de filas n de la matriz A
    NROW = 1:n; % Inicializar el vector de índices de fila
    S = zeros(1, n); % Inicializar el vector de escalas

    % Calcular las escalas S para cada fila
    for i = 1:n
        S(i) = max(abs(A(i, 1:n)));
        if S(i) == 0
            error('No existe una solución única: la fila %d de la matriz está compuesta totalmente por ceros.', i);
        end
    end

    % Proceso de eliminación con pivoteo parcial escalado
    for i = 1:n-1
        % Pivoteo parcial escalado
        [~, p] = max(abs(A(NROW(i:n), i)) ./ S(NROW(i:n)));
        p = p + i - 1;

        if A(NROW(p), i) == 0
            error('No existe una solución única: el elemento pivote es cero.');
        end

        % Intercambiar índices de fila si es necesario
        if i ~= p
            NROW([i, p]) = NROW([p, i]);
        end

        % Eliminación Gaussiana
        for j = i+1:n
            m = A(NROW(j), i) / A(NROW(i), i);
            A(NROW(j), :) = A(NROW(j), :) - m * A(NROW(i), :);
        end
    end

    % Verificar si el último elemento pivote es cero
    if A(NROW(n), n) == 0
        error('No existe una solución única: el último elemento pivote es cero.');
    end

    % Sustitución hacia atrás
    x = zeros(n, 1);
    x(n) = A(NROW(n), end) / A(NROW(n), n);
    for i = n-1:-1:1
        x(i) = (A(NROW(i), end) - A(NROW(i), i+1:n) * x(i+1:n)) / A(NROW(i), i);
    end
end