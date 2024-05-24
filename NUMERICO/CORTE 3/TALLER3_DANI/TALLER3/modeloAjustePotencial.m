function c = modeloAjustePotencial(x, y)
    syms h 
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

    % Define la función potencial utilizando la variable simbólica 'h'
    c = alfa2*h^beta2;

end