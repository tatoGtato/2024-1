function x = partialPivotingElimination(A)
    % Toma como entrada una matriz aumentada A que representa un sistema de
    % ecuaciones lineales y devuelve la solución del sistema.

    [n, m] = size(A);

    if m ~= n + 1
        error('La matriz no tiene las dimensiones correctas.');
    end

    NROW = 1:n; % Inicializar apuntadores de fila

    for i = 1:n-1
        % Pivoteo parcial
        [~, p] = max(abs(A(NROW(i:n), i)));
        p = p + i - 1;

        if A(NROW(p), i) == 0
            error('No existe una solución única.');
        end

        % Intercambiar filas si es necesario
        if NROW(i) ~= NROW(p)
            NCOPY = NROW(i);
            NROW(i) = NROW(p);
            NROW(p) = NCOPY;
        end

        % Eliminación
        for j = i+1:n
            m = A(NROW(j), i) / A(NROW(i), i);
            A(NROW(j), :) = A(NROW(j), :) - m * A(NROW(i), :);
        end
    end

    % Comprobación final
    if A(NROW(n), n) == 0
        error('No existe una solución única.');
    end

    % Sustitución hacia atrás
    x = zeros(n, 1);
    x(n) = A(NROW(n), end) / A(NROW(n), n);

    for i = n-1:-1:1
        x(i) = (A(NROW(i), end) - A(NROW(i), i+1:n) * x(i+1:n)) / A(NROW(i), i);
    end
end