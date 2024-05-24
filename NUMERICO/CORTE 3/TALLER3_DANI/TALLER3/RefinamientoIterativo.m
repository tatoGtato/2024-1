% El método de refinamiento iterativo es una técnica utilizada para mejorar la 
% precisión de una solución aproximada encontrada por un método numérico iterativo, 
% como el método de Jacobi o el método de Gauss-Seidel. Este método se utiliza comúnmente 
% en sistemas de ecuaciones lineales o en la resolución de ecuaciones no lineales.
% 
% En sistemas de ecuaciones lineales, el refinamiento iterativo puede ayudar a mejorar 
% la precisión de las soluciones aproximadas encontradas por métodos como el método de Jacobi, 
% el método de Gauss-Seidel o el método de eliminación gaussiana.
% En la resolución de ecuaciones no lineales, el refinamiento iterativo puede utilizarse para 
% mejorar la precisión de las raíces encontradas por métodos como el método de Newton-Raphson.


function [x, K] = RefinamientoIterativo(A, b, x0, tol)
    iter = 0;

    % Extraemos las matrices diagonal, triangular inferior y triangular superior de A
    D = diag(diag(A));  % Matriz diagonal de A
    L = tril(A, -1);    % Parte triangular inferior de A (sin la diagonal)
    U = triu(A, 1);     % Parte triangular superior de A (sin la diagonal)

    % Calculamos la matriz de iteración T y el vector c usando la matriz diagonal D
    T = inv(D) * (L + U);  % T = D^(-1) * (L + U)
    c = inv(D) * b;        % c = D^(-1) * b

    % Calculamos la primera estimación de la solución
    x = T * x0 + c;

    % Iteramos hasta que la norma relativa de la diferencia sea menor que la tolerancia
    while norm(x - x0, inf) / norm(x, inf) > tol
        % Actualizamos la estimación anterior
        x0 = x;

        % Calculamos el residuo r
        r = b - A * x0;
        
        % Solucionamos el sistema Ay = r
         y = A \ r;  

        % Actualizamos la solución
        x = x0 + y;

        % Contamos las iteraciones
        iter = iter + 1;
    end

    % Imprimimos el número de iteraciones para alcanzar la convergencia
    fprintf('Número de iteraciones para convergencia: %d\n', iter);

    % Calculamos el número de condición de la matriz A
    K = norm(A) * norm(inv(A));
end

