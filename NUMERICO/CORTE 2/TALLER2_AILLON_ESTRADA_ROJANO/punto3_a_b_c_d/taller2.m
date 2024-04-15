%% Punto 3a
h = 0.25;
k = 0.06;
a = 0; 
b = 57; % Prueba y error
alpha = 3; % Condicion inicial
N = (b-a)/h;
TW = RK4(@(t, y) ODE3a(t, y, k), a, b, N, alpha);
plot(TW(:,1),TW(:,2));
title('Drenado de tanque con RK4')
xlabel('t')
ylabel('y')
    
%% Punto 3b y(0) = 1
close all
a = 0;
b = 1;
x = linspace(a,b,1000);
alpha = 1; % condicion inicial
h = 0.05;
N = round((b-a)/h);

% Soluciones analiticas
y1 = @(x) (((x + 2*x.^2 - 2)/2).^2);
y2 = @(x) (((x + 2*x.^2 + 2)/2).^2);

% Euler
TY1 = metodo_euler(@(x,y) ODE3b(x,y), a, b,alpha,N);

% Punto Medio o RK2
TW1 = RK2(@(x,y) ODE3b(x,y),a,b,N,alpha);

% Euler modificado
TY2 = euler_mod(@(x,y) ODE3b(x,y),a,b,alpha,N);

% RK4
TW2 = RK4(@(x,y) ODE3b(x,y),a,b,N,alpha);

% Create plots
t = tiledlayout(3,2);
nexttile
plot(TY1(:,1),TY1(:,2),'--p');
title('Euler');
nexttile
plot(TW1(:,1),TW1(:,2),'--b');
title('Punto Medio');
nexttile
plot(TY2(:,1),TY2(:,2),'--r');
title('Euler Modificado');
nexttile
plot(TW2(:,1),TW2(:,2),'--g');
title('RK4');
nexttile
plot(x,y1(x));
title('Solucion análitica c = -2');
nexttile
plot(x,y2(x));
title('Solución análitica c = 2')
hold off

%% Punto 3c
% Define los parámetros del sistema
m = 20; % Masa
c = 5; % Coeficiente de amortiguamiento
k = 20; % Constante del resorte
h = 0.05; % tamaño de paso

% Condiciones iniciales
x0 = 1; % Posición inicial
v0 = 0; % Velocidad inicial
alpha = [x0; v0]; % Vector de condiciones iniciales

% Intervalo de tiempo
a = 0;
b = 15;
N = 100; % Número de pasos

% Resolver el sistema usando RK4
TW = RK4_system(@(t, Y) sistema(t, Y, m, c, k), a, b, N, alpha);

% Extraer los resultados
T = TW(:,1);
X = TW(:,2);
V = TW(:,3);

for c = [5 40 200]
    TW = RK4_system(@(t, Y) sistema(t, Y, m, c, k), a, b, N, alpha);
    T = TW(:,1);
    X = TW(:,2);
    V = TW(:,3);
    plot(X,T);
    hold on
end


function dYdt = sistema(t,Y,m,c,k)
x = Y(1);
v = Y(2);

% Derivadas
dxdt = v;
dvdt = (-c/m) * v - (k/m) * x;

dYdt = [dxdt; dvdt];
end

%% 3d Euler RK4 Lagrange
clear; 
t = [1950, 1960, 1970, 1980, 1990, 2000];
p = [2555, 3040, 3708, 4454, 5276, 6079];
k = 0.026;
N = 100;
pMax = 12000;

c = lagrangeCoefficients(t,p);
poly = @(t) polyval(c,t);
t_values = linspace(min(t),max(t),1000);
p_values = poly(t_values);


TP_all = [];
TP1_all = [];
for i = 1:length(t)-1
    a = t(i);
    b = t(i+1);
    alpha = p(i);
    
    % Euler
    TP = metodo_euler(@(t,p) k * p * (1 - p/pMax), a, b, alpha, N);
    TP_all = [TP_all; TP];
    % RK4 
    TP1 = RK4(@(t,p) k * p * (1 - p/pMax), a, b, N, alpha);
    TP1_all = [TP1_all; TP1];
end
figure;
hold on;
plot(TP_all(:,1),TP_all(:,2),'Color','red','LineWidth',0.2);
plot(TP1_all(:,1),TP1_all(:,2),'Color','black','LineStyle','--','LineWidth',2);
plot(t_values,p_values,'Color','blue','LineWidth',1);
xlabel('t');
ylabel('P');
title('Euler vs RK4 vs Lagrange');
legend('Euler', 'RK4','Lagrange');
hold off;

