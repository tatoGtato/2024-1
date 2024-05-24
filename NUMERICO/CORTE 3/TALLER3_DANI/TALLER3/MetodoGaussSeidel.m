% En el método de Jacobi, se utilizan los valores x^(k) de la iteración anterior para calcular x(k+1)
% En el método de Gauss-Seidel, se utilizan los valores actualizados x^(k+1)
% tan pronto como están disponibles, es decir, se utilizan los valores recién calculados para 
% actualizar los valores de x durante la misma iteración.

function [x, T, c] = MetodoGaussSeidel(A, b, x0, tol)
    iter = 0;
    % Extraemos las matrices diagonal, triangular inferior y triangular superior de A
    D = diag(diag(A));  % Matriz diagonal
    L = tril(A, -1);    % Parte triangular inferior de A (sin la diagonal)
    U = triu(A, 1);     % Parte triangular superior de A (sin la diagonal)

    % Calculamos la matriz de iteración T y el vector c
    % Evitamos el uso de inv para mejorar la eficiencia y estabilidad numérica
    DL = D + L;
    T = -DL \ U;      
    c = DL \ b;       
    % Inicializamos la estimación de la solución
    x = x0;
    % Iteramos hasta que la norma relativa de la diferencia sea menor que la tolerancia
    while norm(T*x + c - x, inf) / norm(x, inf) > tol
        x0 = x;

        % Actualizamos cada componente de x
        for i = 1:length(A)
            x(i) = T(i, :)*x + c(i);
            iter = iter + 1; 
        end
    end
    fprintf('Número de iteraciones para convergencia: %d\n', iter);
end