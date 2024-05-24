function x = ElimiGausiana(A, b)
    % Obtener el tamaño de la matriz A
    [n, ~] = size(A);
    % Inicializar el vector solución x con ceros
    x = zeros(n, 1);
    % Inicializar el contador de iteraciones
    iter = 0;

    % Fase de eliminación hacia adelante
    for i = 1:n
        iter = iter + 1;
        % Verificar si el pivote es cero
        if A(i, i) == 0
            % Buscar una fila posterior con un elemento no cero en la misma columna
            for k = i+1:n
                iter = iter + 1;
                if A(k, i) ~= 0
                    % Intercambiar las filas i y k en la matriz A
                    A([i k], :) = A([k i], :);
                    % Intercambiar los elementos correspondientes en el vector b
                    b([i k]) = b([k i]);                   
                    break;
                end
            end
        end

        % Realizar la eliminación de Gauss para la columna i
        for j = i+1:n
            % Calcular el multiplicador lambda
            lambda = A(j, i) / A(i, i);
            % Actualizar la fila j de la matriz A
            A(j, i:end) = A(j, i:end) - lambda * A(i, i:end);
            % Actualizar el vector b
            b(j) = b(j) - lambda * b(i);
            iter = iter + 1;
        end
    end

    % Fase de sustitución hacia atrás
    for i = n:-1:1
        % Calcular el valor de x(i) utilizando la sustitución hacia atrás
        x(i) = (b(i) - A(i, i+1:end) * x(i+1:end)) / A(i, i);
        iter = iter + 1;
    end
    disp(['Número total de iteraciones en Eliminación Gaussiana: ', num2str(iter)]);
end
