function c = modeloRazonCrecimiento(x, y)
    syms h 

    grado = 1; % Grado del polinomio para la regresión

    % Realiza una regresión polinómica de grado 1 (lineal) utilizando los inversos de 'x' e 'y'
    test = regresionPolinomial (1./x, 1./y, grado);

    % Invierte el orden de los coeficientes del polinomio (probablemente para que estén en orden de grado decreciente)
    test = flip(test);

    % Extrae los coeficientes del polinomio
    a1 = test(1); % Pendiente de la regresión
    a0 = test(2); % Término independiente de la regresión

    % Calcula los parámetros 'alfa3' y 'beta3' del modelo de razón de crecimiento
    alfa3 = 1/a0;
    beta3 = a1*alfa3;

    % Calcula la función de razón de crecimiento utilizando la variable simbólica 'h'
    c = alfa3*h/(beta3+h);

end
% La función devuelve 'c', que es una expresión simbólica que representa la función de razón de crecimiento en términos de 'h'