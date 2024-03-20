%%PLOT
x = linspace(0,2,10);
plot(x, x - 0.2*sin(x) - 0.5);
axes0


%% BISSECCION
clc
tol = 10^(-4);
f = @(x) x - 0.2*sin(x) - 0.5;
No = 20;
a = 0;
b = 1;


tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 0.6154
%ITERACIONES: 13

%% NEWTON RAPHSON
clc

f = @(x) x - 0.2*sin(x) - 0.5;
fp = @(x) 1 - 0.2*cos(x)
p0 = 0;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 0.6155
%ITERACIONES: 3

%% SECANTE

clc
f = @(x) x - 0.2*sin(x) - 0.5;
p0 = 0;
p1 = 1;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 0.6155
%ITERACIONES: 5

%% FALSA POSICION 

clc
f = @(x) x - 0.2*sin(x) - 0.5;
p0 = 0;
p1 = 1;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 0.6155
%ITERACIONES: 5