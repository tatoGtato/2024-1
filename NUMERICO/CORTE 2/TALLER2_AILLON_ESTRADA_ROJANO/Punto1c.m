clc 
clear all

f = @(x) exp( -(x*x/2) )*(1/(sqrt(2*pi)));
a = -0.5;
b = 0.5;

h = 0.1;
inflX = 0;

for x = a : 0.01 : b
    if (richardson_extrapolation_3puntos(f, x, h) == 0)
       inflX = x
    end
end

f(inflX)


%%
function [derivada] = richardson_extrapolation_3puntos(f, x, h)
    %3 puntos
    derivada = (1/(2*h)) * (-3*f(x) + 4*f(x + h) - f(x + 2*h));

    % Extrapolación de Richardson
    f1 = richardson_aux(f, x, h, 1);
    f2 = richardson_aux(f, x, h, 2);

    % Corrección de la derivada utilizando extrapolación de Richardson
    derivada = (4*f1 - f2) / 3;
end

function [f] = richardson_aux(f, x, h, n)
    % Función auxiliar para la extrapolación de Richardson.
    f = (1/(2^n*h)) * (-f(x + h) + f(x - h));
end