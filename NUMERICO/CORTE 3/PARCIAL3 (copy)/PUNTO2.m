clc 
clear all

A = [4  1  1  0  1;
    -1 -3  1  1  0;
     2  1  5 -1 -1;
    -1 -1 -1  4  0;
     0  2 -1  1  4];

b = [6;
     6;
     6;
     6;
     6];

x0 = [0;
      0;
      0;
      0;
      0];

tol = 10^(-3);

[x, T, c] = MetodoJacobi(A, b, x0, tol)

%Luego de 15 iteraciones el resultado del sistema con Jacobi fue
   % x1 = 0.7871
   % x2 = -1.0030
   % x3 = 1.8660
   % x4 = 1.9124
   % x5 = 1.9896

%% FUNCIONES UTILIZADAS

function [x, T, c] = MetodoJacobi(A, b, x0, tol)
    iter = 0; 

    D = diag(diag(A)); % Extraer la diagonal de la matriz A
    L = tril(A, -1); % Obtener la parte triangular inferior de A sin la diagonal
    U = triu(A, 1); % Obtener la parte triangular superior de A sin la diagonal

    % En lugar de invertir D, usamos la división de matrices de MATLAB
    T = -D \ (L + U); % Calcular la matriz T para la iteración
    c = D \ b; % Calcular el vector c para la iteración

    x = T*x0 + c; % Realizar la primera iteración

    while norm(x - x0, inf) / norm(x, inf) > tol % Comprobar la condición de tolerancia
        x0 = x; % Actualizar el vector x0 con el resultado de la iteración anterior
        x = T*x + c; % Realizar la siguiente iteración
        iter = iter + 1; 
    end 
    fprintf('Número de iteraciones para convergencia: %d\n', iter); 
end