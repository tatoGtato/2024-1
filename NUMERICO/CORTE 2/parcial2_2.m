clc
clear all


g = 9.8;
R = 6.37*10^(6);

a = 0;
b = 500;
N = 10;
alpha = 1500;

TW = RK4_system(@(t, Y) sistema(t, Y, R, g), a, b, N, alpha)
T = TW(:,1);
x = TW(:,2);

% Top plot
nexttile
plot(T, x);
title('Aproximacion')
xlim([0 500])




%%
%FUNCIONES UTILIZADAS

function TW = RK4_system(f, a, b, N, alpha, R, g)
    % Paso 1: Inicialización de variables
    h = (b - a) / N;  % Tamaño del paso
    T = a:h:b;  % Vector de tiempos desde a hasta b con N+1 puntos
    Y = zeros(2, N+1);  % Matriz para almacenar las aproximaciones de las variables de estado
    Y(:,1) = alpha;  % Valor inicial de Y

    % Paso 2: Bucle para N iteraciones
    for i = 1:N
        % Paso 3: Calcular K1, K2, K3, K4
        K1 = h * f(T(i), Y(:,i));
        K2 = h * f(T(i) + h/2, Y(:,i) + K1/2);
        K3 = h * f(T(i) + h/2, Y(:,i) + K2/2);
        K4 = h * f(T(i) + h, Y(:,i) + K3);

        % Paso 4: Calcular el nuevo valor de Y
        Y(:,i+1) = Y(:,i) + (K1 + 2*K2 + 2*K3 + K4) / 6;
    end

    % Combinar T y Y en una matriz TW
    TW = [T' Y'];
end

function dYdt = sistema(t, Y, R, g)
    x = Y(1);
    v = Y(2);

    % Derivadas
    dxdt = v;
    dvdt =  -g*((R^2)/(R+x)^2);
    
    dYdt = [dxdt; dvdt];
end

