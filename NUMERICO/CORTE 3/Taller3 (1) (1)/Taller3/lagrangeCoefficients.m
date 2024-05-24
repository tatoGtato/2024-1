function P = lagrangeCoefficients(X, y)
    syms x
    n = length(X); % Número de puntos
    P = 0; % Inicializa el polinomio resultante
    
    % Itera sobre cada punto de datos
    for i = 1:n
        % Inicializa el polinomio base de Lagrange Li(x) = 1
        Li = 1;
        
        % Construye el polinomio base de Lagrange para el punto i
        for j = 1:n
            if j ~= i
                Li = Li * ((x-X(j))/(X(i)-X(j)));
            end
        end
        % Ajusta el tamaño del polinomio Li para que coincida con P
        Li = [zeros(1, length(P) - length(Li)), Li];
        
        % Suma el polinomio base ponderado por yi al polinomio resultante
        P = P + y(i) * Li;
    end
end