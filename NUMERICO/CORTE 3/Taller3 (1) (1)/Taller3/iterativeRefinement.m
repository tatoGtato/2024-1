function [xx, COND] = iterativeRefinement(A, b, N, TOL, t)
    % Paso 3: Resolver el sistema Ax = b por eliminación Gaussiana
    [L, U] = luFactorization(A);
    y = forwardSubstitution(L, b);
    x = backwardSubstitution(U, y);

    % Paso 4: Inicializar k
    k = 1;
    
    % Paso 5: Iterar hasta alcanzar el número máximo de iteraciones
    while k <= N
        % Paso 6: Calcular r_i
        r = b - A * x;

        % Paso 7: Resolver Ay r usando eliminación Gaussiana
        y = forwardSubstitution(L, r);
        y = backwardSubstitution(U, y);

        % Paso 10: Calcular xx_i
        xx = x + y;

        % Paso 13: Calcular COND en la primera iteración
        if k == 1
            COND = norm(y, inf) / norm(x, inf) * 10^t;
        end

        % Paso 16: Comprobar si se alcanzó la tolerancia
        if norm(x - xx, inf) < TOL
            % Paso 17: Salida exitosa
            return;
        end

        % Paso 22: Incrementar k
        k = k + 1;

        % Paso 23: Actualizar x
        x = xx;
    end

    % Paso 27: Salida no exitosa (excedió el número máximo de iteraciones)
    return;
end

function y = forwardSubstitution(L, b)
    % Resuelve el sistema triangular inferior Ly = b
    n = length(b);
    y = zeros(n, 1);
    y(1) = b(1) / L(1, 1);
    for i = 2:n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end
end

function x = backwardSubstitution(U, y)
    % Resuelve el sistema triangular superior Ux = y
    n = length(y);
    x = zeros(n, 1);
    x(n) = y(n) / U(n, n);
    for i = n-1:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end
