% El método de Jacobi es un algoritmo iterativo utilizado para resolver sistemas de ecuaciones lineales.
% Este método es especialmente útil cuando se tiene un sistema de ecuaciones lineales grandes y dispersas, 
% donde la matriz de coeficientes es grande y tiene muchos ceros.

function [x, T, c] = MetodoJacobi(A, b, x0, tol)
    iter = 0; 

    D = diag(diag(A)); % Extraer la diagonal de la matriz A
    L = tril(A, -1); % Obtener la parte triangular inferior de A sin la diagonal
    U = triu(A, 1); % Obtener la parte triangular superior de A sin la diagonal

    % En lugar de invertir D, usamos la división de matrices de MATLAB
    T = -D \ (L + U); % Calcular la matriz T para la iteración
    c = D \ b; % Calcular el vector c para la iteración

    x = T*x0 + c; % Realizar la primera iteración

    while norm(x - x0, inf) / norm(x, inf) > tol % Comprobar la condición de tolerancia
        x0 = x; % Actualizar el vector x0 con el resultado de la iteración anterior
        x = T*x + c; % Realizar la siguiente iteración
        iter = iter + 1; 
    fprintf('Número de iteraciones para convergencia: %d\n', iter); 
end