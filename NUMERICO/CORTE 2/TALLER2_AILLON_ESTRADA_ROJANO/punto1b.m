clear all
clc

f = @(x) cos(x);
x = pi/4;
h = pi/4;

richardson_extrapolation_3puntos(f, x, h)

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