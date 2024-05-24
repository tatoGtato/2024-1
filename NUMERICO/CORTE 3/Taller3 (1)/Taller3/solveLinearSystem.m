function [x, uniqueSol] = solveLinearSystem(A)
    [n, m] = size(A);
    
    if m ~= n + 1
        error('La matriz proporcionada no es una matriz aumentada válida para un sistema n x n.');
    end
    
    % Paso 1: Proceso de eliminación
    for i = 1:n-1
        % Paso 2: Buscar p tal que api != 0
        p = find(A(i:n, i) ~= 0, 1) + i - 1;
        
        if isempty(p)
            x = [];
            uniqueSol = false;
            disp('No existe una solución única.');
            return;
        end
        
        % Paso 3: Intercambiar filas si es necesario
        if p ~= i
            A([i p], :) = A([p i], :);
        end
        
        % Paso 4: Eliminación para las filas j = i+1 a n
        for j = i+1:n
            mji = A(j, i) / A(i, i);
            A(j, :) = A(j, :) - mji * A(i, :);
        end
    end
    
    % Paso 7: Verificar si la solución es única
    if A(n, n) == 0
        x = [];
        uniqueSol = false;
        disp('No existe una solución única.');
        return;
    end
    
    % Paso 8: Sustitución hacia atrás
    x = zeros(n, 1);
    x(n) = A(n, m) / A(n, n);
    
    % Paso 9: Continuar con la sustitución hacia atrás
    for i = n-1:-1:1
        x(i) = (A(i, m) - A(i, i+1:n) * x(i+1:n)) / A(i, i);
    end
    
    uniqueSol = true;
end