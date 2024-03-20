%%PLOT

x = linspace(0,2,10);
plot(x, -2*x.^6 - 1.5*x.^4 +10*x+2);
axes0

%% a) FALSA POSICION 
clc

f = @(x) -2*x^6 - 1.5*x^4 +10*x+2;
fp = @(x) -12*x^5 - 6*x^3 + 10;
tol = 10^(-4);
No = 20;

p0 = 0;
p1 = 1;

[p,i] = falsaPosicion(fp, p0, p1, tol, No)

max = f(p)

%% b) NEWRON-RAPHSON
clc

f = @(x) -2*x^6 - 1.5*x^4 +10*x+2;
fp = @(x) -12*x^5 - 6*x^3 + 10;
fp2 = @(x) -60*x^4 - 18*x^2;
tol = 10^(-4);
No = 20;

p0 = 1;

[p, i] = newtonRaphson(fp, fp2, p0, tol, No)

max = f(p)

%% c) SECANTE
clc

f = @(x) -2*x^6 - 1.5*x^4 +10*x+2;
fp = @(x) -12*x^5 - 6*x^3 + 10;
tol = 10^(-4);
No = 20;

p0 = 0;
p1 = 1;

[p,i] = secante(fp, p0, p1, tol, No)

max = f(p)