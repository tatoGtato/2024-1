%%PLOT
x = linspace(-10,10);
plot(x, exp(x).*cos(x) - x.^2 + 3*x);
ylim([-5 5])
axes0

%% BISSECCION [-1, 8]
clc
tol = 10^(-4);
f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
No = 20;

a = -1;
b = 8;
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 4.7838
%ITERACIONES: 16

%% BISSECCION [1, 2]
clc
tol = 10^(-4);
f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
No = 20;

a = 1;
b = 2;
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 4.7838
%ITERACIONES: 16

%% BISSECCION [-1, 0]
clc
tol = 10^(-4);
f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
No = 20;

a = -1;
b = 0;
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 4.7838
%ITERACIONES: 16

%% NEWTON RAPHSON P0 = -2
clc

f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
fp = @(x) exp(x).*cos(x) - exp(x).*sin(x) - 2*x + 3;
p0 = -2;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 1.073
%ITERACIONES: 5

%% NEWTON RAPHSON P0 = 2
clc

f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
fp = @(x) exp(x).*cos(x) - exp(x).*sin(x) - 2*x + 3;
p0 = 2;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 1.073
%ITERACIONES: 5

%% NEWTON RAPHSON P0 = 5
clc

f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
fp = @(x) exp(x).*cos(x) - exp(x).*sin(x) - 2*x + 3;
p0 = 5;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 1.073
%ITERACIONES: 5

%% SECANTE P0 = 5, p1 = 6

clc
f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
p0 = 5;
p1 = 6;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 1.073
%ITERACIONES: 7

%% FALSA POSICION 

clc
f = @(x) exp(x).*cos(x) - x.^2 + 3*x;
p0 = 5;
p1 = 6;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 1.037
%ITERACIONES: 10


