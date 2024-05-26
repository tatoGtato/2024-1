clc
clear all
syms X

xi = [4.0,  4.2,  4.5,  4.7,  5.1,  5.5,  5.9,  6.3,  6.8,  7.1];
yi = [102.56,  113.18,  130.11,  142.05,  167.53,  195.14,  224.87,  256.73,  299.50,  326.72];

%% MINIMOS CUADRADOS GRADO 1
c = regresionPolinomial(xi, yi, 1);
test = flip(c);
aprox1 = vpa(poly2sym(test, X))

error_estandar = sqrt(sum((subs(aprox1, X, xi) - yi).^2) / (length(yi) - 2))

%POLINOMIO: 72.0845*X - 194.1382

%% MINIMOS CUADRADOS GRADO 2
c = regresionPolinomial(xi, yi, 2);
test = flip(c);
aprox2 = vpa(poly2sym(test, X))

error_estandar = sqrt(sum((subs(aprox2, X, xi) - yi).^2) / (length(yi) - 2))

%POLINOMIO: 6.6182*X^2 - 1.1435*X + 1.2356

%% DE LA FORMA b*e^ax
[b,a] = regresion_exponencial(xi,yi);
aproxEXP = b*exp(a*X)

error_estandar = sqrt(sum((subs(aproxEXP, X, xi) - yi).^2) / (length(yi) - 2))

%POLINOMIO: 24.2588*exp(0.3724*X)

%% DE LA FORMA b*x^a
[b,a] = regresion_potencial(xi,yi);
aproxPOT = b*X^a

error_estandar = sqrt(sum((subs(aproxPOT, X, xi) - yi).^2) / (length(yi) - 2))

%POLINOMIO: 6.2390*X^2.0195

%% GRAFICAS

hold on

plot(xi, yi, "o")
plot(xi, subs(aprox1, X, xi))
plot(xi, subs(aprox2, X, xi))
plot(xi, subs(aproxEXP, X, xi))
plot(xi, subs(aproxPOT, X, xi))

legend('Datos','Minimos cuadrados de grado 1','Minimos cuadrados de grado 2', 'Exponencial', 'Potencial')

hold off

%% FUNCIONES UTILIZADAS

function [a0,a1] = regresion_lineal(xi,yi)
    n = length(xi);
    a1 = (length(xi)*sum(xi.*yi)-sum(xi)*sum(yi))/(n*sum(xi.^2)-(sum(xi))^2);
    a0 = mean(yi)-a1*mean(xi);
end

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

function [beta2,alfa2] = modeloAjustePotencial(x, y)
    grado = 1; 
    % Aplica regresión polinomial a los logaritmos de base 10 de 'x' e 'y'
    % Esto se hace porque si y = alpha * x^beta, entonces log10(y) = log10(alpha) + beta * log10(x)
    % lo cual es una relación lineal entre log10(y) y log10(x)
    test = regresionPolinomial (log10(x), log10(y), grado);

    % Invierte el orden de los coeficientes del polinomio para que estén en orden de grado decreciente
    test = flip(test);

    % Extrae los coeficientes del polinomio
    a1 = test(1); % Pendiente de la regresión, que corresponde a 'beta' en el modelo potencial
    a0 = test(2); % Término independiente de la regresión, que corresponde a 'log10(alpha)'

    % Calcula los parámetros 'alfa2' y 'beta2' del modelo potencial
    beta2 = a1; % 'beta' en el modelo potencial
    alfa2 = 10^(a0); % 'alpha' en el modelo potencial, calculado como 10 elevado al término independiente de la regresión

end

function [a0,a1] = regresion_exponencial(xi,yi)
    x2 = xi;
    y2 = log(yi);
    
    [A0,A1] = regresion_lineal(x2,y2);

    a1 = A1;
    a0 = exp(A0);
    %funct(x)= A0*exp(A1*x);
end

%está estructurado como una función que toma una matriz matriz y un vector z, y devuelve el vector solución x.
function x = eliminacionGaussianaPivoteoEscalado(matriz, z)
    % Obtiene el tamaño de la matriz de coeficientes y el vector de resultados
    tamano = size(matriz);
    tamanoZ = size(z);
    
    % Calcula los factores de escalamiento como el máximo valor absoluto por fila
    scalingFactors = max(abs(matriz), [], 2); 
    
    % Proceso de eliminación gaussiana con pivoteo escalado
    for i = 1:tamano(2)
        % Encuentra el valor máximo escalado en la columna actual y su fila pivote
        [maxValue, pivotRow] = max(abs(matriz(i:end, i)) ./ scalingFactors(i:end));
        pivotRow = pivotRow + i - 1;  % Ajusta el índice de la fila pivote
        
        % Intercambia la fila actual con la fila pivote si es necesario
        if pivotRow ~= i
            matriz([i, pivotRow], :) = matriz([pivotRow, i], :);
            z([i, pivotRow]) = z([pivotRow, i]);
        end
        
        % Realiza la eliminación para hacer ceros debajo de la diagonal principal
        for j = (i + 1):tamano(1)
            factor = matriz(j, i) / matriz(i, i);  % Calcula el factor de eliminación
            matriz(j, :) = matriz(j, :) - factor * matriz(i, :);  % Actualiza la fila j
            z(j) = z(j) - factor * z(i);  % Actualiza el vector de resultados
        end
    end
    
    % Inicializa el vector de solución x con ceros
    x = zeros(1, tamanoZ(1));
    
    % Proceso de sustitución hacia atrás para obtener la solución del sistema
    for i = tamanoZ(1):-1:1
        x(i) = z(i);  % Inicializa x(i) con el valor correspondiente en z
        for j = (i + 1):tamanoZ(1)
            x(i) = x(i) - matriz(i, j) * x(j);  % Resta los términos conocidos
        end
        x(i) = x(i) / matriz(i, i);  % Divide por el coeficiente diagonal para obtener x(i)
    end
end

function [a0,a1] = regresion_potencial(xi,yi)
    x2 = log10(xi);
    y2 = log10(yi);
    
    [A0,A1] = regresion_lineal(x2,y2);
    a1 = A1;
    a0 = 10^A0;
    %funct(x)= A0*x^A1;
end