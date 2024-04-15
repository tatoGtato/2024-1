function TY = metodo_euler(f, a, b, alpha, N)
    % Inicialización de variables
    h = (b - a) / N; % Tamaño del paso
    t = a:h:b; % Vector de tiempo desde a hasta b con incrementos de h
    y = zeros(1, N+1); % Vector de soluciones inicializado a cero
    y(1) = alpha; % Condición inicial

    % Método de Euler
    for i = 1:N
        y(i+1) = y(i) + h * f(t(i), y(i)); % Actualización de y
    end
    TY = [t' y'];
end
