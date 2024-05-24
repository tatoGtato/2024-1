function [L, U] = luFactorization(A)
    n = size(A, 1); % Número de filas/columnas de A
    L = eye(n); % Inicializa L como la matriz identidad
    U = zeros(n); % Inicializa U como una matriz de ceros

    % Paso 3
    if A(1,1) == 0
        error('Factorización imposible');
    else
        L(1,1) = 1; % Elegimos l11 como 1 para simplificar
        U(1,1) = A(1,1); % u11 = a11
    end

    % Primer bucle: Establece la primera fila de U y la primera columna de L
    for j = 2:n
        U(1,j) = A(1,j) / L(1,1);
        L(j,1) = A(j,1) / U(1,1);
    end

    % Segundo bucle: Cálculo de los elementos restantes de L y U
    for i = 2:n-1
        for k = 1:i-1
            L(i,i) = 1; % Asumimos lii como 1
            U(i,i) = A(i,i) - sum(L(i,1:k) .* U(1:k,i)');
        end

        if U(i,i) == 0
            error('Factorización imposible');
        end

        for j = i+1:n
            U(i,j) = (A(i,j) - sum(L(i,1:i-1) .* U(1:i-1,j)')) / L(i,i);
            L(j,i) = (A(j,i) - sum(L(j,1:i-1) .* U(1:i-1,i)')) / U(i,i);
        end
    end

    % Último elemento de U (y L ya está completo con 1s en la diagonal)
    L(n,n) = 1;
    U(n,n) = A(n,n) - sum(L(n,1:n-1) .* U(1:n-1,n)');

    if U(n,n) == 0
        error('La matriz A es singular, no existe una solución única');
    end
end