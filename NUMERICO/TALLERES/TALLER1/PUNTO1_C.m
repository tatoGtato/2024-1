%%PLOT
x = linspace(-4,10,10);
plot(x, exp(x/2) - x.^2 - 3*x);
axes0


%% BISSECCION [-4, 10]
clc
tol = 10^(-4);
f = @(x) exp(x/2) - x^2 - 3*x
No = 20;
a = -4;
b = 10;


tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 9.5856
%ITERACIONES: 17

%% BISSECCION [0, 2]
clc
tol = 10^(-4);
f = @(x) exp(x/2) - x^2 - 3*x
No = 20;
a = 0;
b = 2;


tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: 9.5856
%ITERACIONES: 14


%% BISSECCION [-4, 1]
clc
tol = 10^(-4);
f = @(x) exp(x/2) - x^2 - 3*x
No = 20;
a = -4;
b = 1;


tic
[p, i] = biseccion(f, a, b, tol, No)
toc

%RESPUESTA: -3.0702
%ITERACIONES: 15

%% NEWTON RAPHSON p0 = -4
clc

f = @(x) exp(x/2) - x^2 - 3*x;
fp = @(x) (exp(x/2)/2) - 2*x - 3
p0 = -4;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: -3.0702
%ITERACIONES: 4

%% NEWTON RAPHSON p0 = 8
clc

f = @(x) exp(x/2) - x^2 - 3*x;
fp = @(x) (exp(x/2)/2) - 2*x - 3
p0 = 8;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 9.5856
%ITERACIONES: 7

%% NEWTON RAPHSON p0 = 0
clc

f = @(x) exp(x/2) - x^2 - 3*x;
fp = @(x) (exp(x/2)/2) - 2*x - 3
p0 = 0;
tol = 10^(-4);
No = 20;

tic
[p, i] = newtonRaphson(f, fp, p0, tol, No)
toc 

%RESPUESTA: 0.3560
%ITERACIONES: 4

%% SECANTE p0 = -4, p1 = 10

clc
f = @(x) exp(x/2) - x^2 - 3*x;
p0 = -4;
p1 = 10;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: -3.0702
%ITERACIONES: 8

%% SECANTE p0 = 8, p1 = 10

clc
f = @(x) exp(x/2) - x^2 - 3*x;
p0 = 8;
p1 = 10;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 9.5856
%ITERACIONES: 6


%% SECANTE p0 = 0, p1 = 2

clc
f = @(x) exp(x/2) - x^2 - 3*x;
p0 = 0;
p1 = 2;
tol = 10^(-4);
No = 20;

tic
[p,i] = secante(f, p0, p1, tol, No)
toc

%RESPUESTA: 0.3560
%ITERACIONES: 7

%% FALSA POSICION p0 = -4, p1 = 10

clc
f = @(x) exp(x/2) - x^2 - 3*x;
p0 = -4;
p1 = 10;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: -3.0702
%ITERACIONES: 10

%% FALSA POSICION p0 = 8, p1 = 10

clc
f = @(x) exp(x/2) - x^2 - 3*x;
p0 = 8;
p1 = 10;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 9.5856
%ITERACIONES: 8

%% FALSA POSICION p0 = 0, p1 = 2

clc
f = @(x) exp(x/2) - x^2 - 3*x;
p0 = 0;
p1 = 2;
tol = 10^(-4);
No = 20;

tic
[p,i] = falsaPosicion(f, p0, p1, tol, No)
toc

%RESPUESTA: 0.3560
%ITERACIONES: 10