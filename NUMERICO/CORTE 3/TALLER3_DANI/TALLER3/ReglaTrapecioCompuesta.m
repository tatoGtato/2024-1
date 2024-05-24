function integral = ReglaTrapecioCompuesta(f, a, b, n)
    h = (b - a) / n;
    xj = a + (1:(n-1)) * h;

    integral = (h / 2) * (f(a) + 2 * sum(f(xj)) + f(b));
end
% f: la función a integrar, que debe ser una función handle en MATLAB.
% 
% lower_lim: el límite inferior del intervalo de integración.
% 
% upper_lim: el límite superior del intervalo de integración.
% 
% n o num_subintervals: el número de subintervalos en los que dividir el intervalo de integración, 
% o el número de puntos de evaluación incluyendo los extremos. En el código proporcionado, 
% se utiliza un valor fijo de 1000, lo que sugiere que se está especificando el número de puntos de evaluación.