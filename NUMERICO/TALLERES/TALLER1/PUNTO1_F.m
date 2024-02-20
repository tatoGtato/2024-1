%%PLOT
x = linspace(-4,6);
plot(x,exp(x) - 4*x.^2 - 8*x);
axes0


%% BISSECCION [-4,6]
clc
tol = 10^(-4);
f = @(x) exp(x) - 4*x^2 - 8*x;
No = 20;

a = -4;
b = 6;
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 4.9108
%ITERACIONES: 16

%% BISSECCION [0,1]
clc
tol = 10^(-4);
f = @(x) exp(x) - 4*x^2 - 8*x;
No = 20;

a = 0;
b = 1;
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 0.1340
%ITERACIONES: 14


%% BISSECCION [-3,-1]
clc
tol = 10^(-4);
f = @(x) exp(x) - 4*x^2 - 8*x;
No = 20;

a = -3;
b = -1;
fa = f(a);
i = 0;

tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: -2.0165
%ITERACIONES: 14

%% NEWTON RAPHSON p0 = -4
clc

f = @(x) exp(x) - 4*x^2 - 8*x;
fp = @(x) exp(x) - 8*x - 8;
p0 = -4;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: -2.0165
%ITERACIONES: 5

%% NEWTON RAPHSON p0 = 0
clc

f = @(x) exp(x) - 4*x^2 - 8*x;
fp = @(x) exp(x) - 8*x - 8;
p0 = 0;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: -2.0165
%ITERACIONES: 8


%% NEWTON RAPHSON p0 = 4
clc

f = @(x) exp(x) - 4*x^2 - 8*x;
fp = @(x) exp(x) - 8*x - 8;
p0 = 4;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 4.9108
%ITERACIONES: 8

%% SECANTE p = -4, p1 = 6

clc
f = @(x) exp(x) - 4*x^2 - 8*x;
p0 = -4;
p1 = 6;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: -2.0165
%ITERACIONES: 8

%% SECANTE p0 = 0, p1 = 1

clc
f = @(x) exp(x) - 4*x^2 - 8*x;
p0 = 0;
p1 = 1;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 0.1339
%ITERACIONES: 6

%% SECANTE p0 = 4, p1 = 6

clc
f = @(x) exp(x) - 4*x^2 - 8*x;
p0 = 4;
p1 = 6;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 4.9108
%ITERACIONES: 9

%% FALSA POSICION p0 = -4, p1 = 6

clc
f = @(x) exp(x) - 4*x^2 - 8*x;
p0 = -4;
p1 = 0;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: -2.0164
%ITERACIONES: 19

%% FALSA POSICION  p0 = 0, p1 = 1

clc
f = @(x) exp(x) - 4*x^2 - 8*x;
p0 = 0;
p1 = 1;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 0.1339
%ITERACIONES: 10

%% FALSA POSICION p0 = 4, p1 = 6

clc
f = @(x) exp(x) - 4*x^2 - 8*x;
p0 = 4;
p1 = 6;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 4.9107
%ITERACIONES: 17
