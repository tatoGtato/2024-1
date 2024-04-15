function TW = RK2(f, a, b, N, alpha)
    % Paso 1: Inicialización de variables
    h = (b - a) / N;  % Tamaño del paso
    T = a:h:b;  % Vector de tiempos desde a hasta b con N+1 puntos
    W = zeros(1, N+1);  % Vector para almacenar las aproximaciones de w
    W(1) = alpha;  % Valor inicial de w (aproximación de y)

    % Paso 2: Bucle para N iteraciones
    for i = 1:N
        % Paso 3: Calcular K1, K2, K3, K4
        K1 = h * f(T(i), W(i));
        K2 = h * f(T(i) + h/2, W(i) + K1/2);

        % Paso 4: Calcular el nuevo valor de w
        W(i+1) = W(i) + (K1 + K2) / 2;
    end

    % Combinar T y W en una matriz TW
    TW = [T' W'];
end