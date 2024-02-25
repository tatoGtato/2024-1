%%PLOT
x = linspace(-1,3);
plot(x,0.53*x.^3 + x.^2 - 2*x - 5);
axes0



%% BISSECCION
clc
tol = 10^(-4);
f = @(x) 0.53*x^3 + x^2 - 2*x - 5;
No = 20;

a = 1;
b = 3;1,
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 2.0871
%ITERACIONES: 14

%% NEWTON RAPHSON
clc

f = @(x) 0.53*x^3 + x^2 - 2*x - 5;
fp = @(x) 1.59*x^2 + 2*x - 2;
p0 = 1;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 2.0871
%ITERACIONES: 7

%% SECANTE

clc
f = @(x) 0.53*x^3 + x^2 - 2*x - 5;
p0 = 1;
p1 = 3;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 2.0871
%ITERACIONES: 8

%% FALSA POSICION 

clc
f = @(x) 0.53*x^3 + x^2 - 2*x - 5;
p0 = 1;
p1 = 3;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 2.0871
%ITERACIONES: 11
