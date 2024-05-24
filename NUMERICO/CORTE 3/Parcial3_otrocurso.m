% Año
syms X

an = [1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020];

% Número de Casos Diagnosticados de diabetes tipo 2 (en millones)
di = [30, 50, 70, 100, 150, 200, 250, 350, 450, 550, 700];

hold on

plot(an, di, 'o');

[a,b] = regresion_exponencial(an,di);

reg = a*exp(b*X);

x1 = linspace(an(1), an(length(anio)), 100);
y1 = subs(reg, X, x1)
plot(x1,y1, '-');

data = table(an', subs(reg, X, an)', 'VariableNames', {'Años', 'polinomio_values'});

Data_1980 = [];
mes = [];

%TABLA 1980
for i = 1 : 12
    mes_ind = i/12
    dat = subs(reg, X, 1980+mes_ind);
    Data_1980(i) = dat; 
    mes(i) = 1980+mes_ind
end

data = table(mes', Data_1980', 'VariableNames', {'Meses de 1980', 'polinomio_values'});

plot(mes, Data_1980, '+');

hold off

%%
% Año
clc 
clear all
syms X


an = [1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020];

% Número de Casos Diagnosticados de diabetes tipo 2 (en millones)
di = [14.5, 15.0, 15.3, 18.4, 23.3, 25.1, 30.5, 34.3, 35.7, 39.6, 42.4];

hold on

plot(an, di, 'o');

[m,b] = ajustePotencias(an, di);

Y_ajuste_saturacion = m*X^b;

x1 = linspace(an(1), an(length(an)), 100);
y1 = subs(Y_ajuste_saturacion, X, x1)
plot(x1,y1, '-');

data = table(an', subs(Y_ajuste_saturacion, X, an)', 'VariableNames', {'Años', 'polinomio_values'});

Data_1980 = [];
mes = [];

%TABLA 1980
for i = 1 : 12
    mes_ind = i/12
    dat = subs(Y_ajuste_saturacion, X, 1980+mes_ind);
    Data_1980(i) = dat; 
    mes(i) = 1980+mes_ind
end

data = table(mes', Data_1980', 'VariableNames', {'Meses de 1980', 'polinomio_values'});

plot(mes, Data_1980, '+');

hold off

%%
function [a0,a1] = regresion_exponencial(xi,yi)
    x2 = xi;
    y2 = log(yi);
    
    [m,b] = ajusteLinealMinimosCuadrados(x2,y2);

    a1 = m;
    a0 = exp(b);
    %funct(x)= A0*exp(A1*x);
end
% 

function [m, b] = ajusteLinealMinimosCuadrados(X, Y)
    % Número de puntos de datos
    n = length(X);

    % Calcular las sumatorias necesarias
    sum_xy = sum(X .* Y);
    sum_x = sum(X);
    sum_y = sum(Y);
    sum_x_squared = sum(X.^2);

    % Calcular la pendiente (m) y la ordenada al origen (b)
    m = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x^2);
    b = (sum_y - m * sum_x) / n;
end

function [a, b] = ajusteCrecimientoSaturacion(X, Y)
    modelo = @(p, x) p(1) .* x ./ (p(2) + x);

    % Estimación inicial de los parámetros
    p0 = [1, 1];

    % Ajustar los parámetros utilizando mínimos cuadrados
    parametros = lsqcurvefit(modelo, p0, X, Y);

    % Extraer los valores ajustados de los parámetros
    a = parametros(1);
    b = parametros(2);
end

function [a0,a1] = regresion_potencial(xi,yi)
    x2 = log10(xi);
    y2 = log10(yi);
    
    [A0,A1] = regresion_lineal(x2,y2);
    a1 = A1;
    a0 = 10^A0;
    %funct(x)= A0*x^A1;
end

function [a0,a1] = regresion_lineal(xi,yi)
    n = length(xi);
    a1 = (length(xi)*sum(xi.*yi)-sum(xi)*sum(yi))/(n*sum(xi.^2)-(sum(xi))^2);
    a0 = mean(yi)-a1*mean(xi);
end

