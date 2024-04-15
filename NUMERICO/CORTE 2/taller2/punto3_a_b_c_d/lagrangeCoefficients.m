function P = lagrangeCoefficients(x, y)
    n = length(x); % Número de puntos
    P = zeros(1, n); % Inicializa el polinomio resultante
    
    % Itera sobre cada punto de datos
    for i = 1:n
        % Inicializa el polinomio base de Lagrange Li(x) = 1
        Li = 1;
        
        % Construye el polinomio base de Lagrange para el punto i
        for j = 1:n
            if j ~= i
                % Crea un polinomio (x - xj) para el término actual
                p = [1, -x(j)];
                % Divide por (xi - xj)
                Li = conv(Li, p) / (x(i) - x(j));
            end
        end
        
        % Ajusta el tamaño del polinomio Li para que coincida con P
        Li = [zeros(1, length(P) - length(Li)), Li];
        
        % Suma el polinomio base ponderado por yi al polinomio resultante
        P = P + y(i) * Li;
    end
end
