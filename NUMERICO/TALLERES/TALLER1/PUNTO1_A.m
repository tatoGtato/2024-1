%% BISSECCION
clc
tol = 10^(-4);
f = @(x) exp(x) - 4 + x;
No = 20;

a = 0;
b = 2;
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 1.073
%ITERACIONES: 14

%% NEWTON RAPHSON
clc

f = @(x) exp(x) - 4 + x;
fp = @(x) exp(x) + 1;
p0 = 0;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 1.073
%ITERACIONES: 5

%% SECANTE

clc
f = @(x) exp(x) - 4 + x;
p0 = 0;
p1 = 2;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 1.073
%ITERACIONES: 7

%% FALSA POSICION 

clc
f = @(x) exp(x) - 4 + x;
p0 = 0;
p1 = 2;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 1.037
%ITERACIONES: 10





