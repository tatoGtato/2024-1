function c = regresionPolinomial(x, y, grado)
    % x: vector de valores de la variable independiente
    % y: vector de valores de la variable dependiente
    % grado: grado del polinomio que se quiere ajustar
    % c: vector de coeficientes del polinomio ajustado

    % Se obtiene la cantidad de datos
    n = length(x);
    
    % Inicializar el vector de sumas de potencias de x
    potenciasDeXi = zeros(1, grado + 1);
    
    % Calcular las sumas de las potencias de x
    for i = 1:grado + 1
        for j = 1:n
            potenciasDeXi(i) = potenciasDeXi(i) + x(j)^(i - 1);
        end
    end
    
    % Inicializar la matriz del sistema de ecuaciones lineales
    matriz = zeros(grado + 1, grado + 1);
    % Inicializar el vector del lado derecho del sistema de ecuaciones
    b = zeros(1, grado + 1);
    
    % Llenar la matriz y el vector b con las sumas correspondientes
    for i = 1:grado + 1
        for j = 1:grado + 1
            for k = 1:n
                matriz(i, j) = matriz(i, j) + x(k)^(i + j - 2);
            end
        end
        
        for j = 1:n
            b(i) = b(i) + x(j)^(i - 1) * y(j);
        end
    end
    
    % Resolver el sistema de ecuaciones lineales usando eliminación gaussiana con pivoteo escalado
    c = eliminacionGaussianaPivoteoEscalado(matriz, b');
    % Alternativamente, se podría usar la división matricial de MATLAB para resolver el sistema
    % c = matriz \ b';
end

% Se utiliza un método de ajuste de curva, como el método de mínimos cuadrados, 
% para encontrar los coeficientes del polinomio que mejor se ajustan a los datos observados. 
% El objetivo es minimizar la suma de los cuadrados de las diferencias entre los valores 
% predichos por el polinomio y los valores observados.