%%
%ANALITICO
clear all
clc

f = @(x) (6 + 3*cos(x));
x0 = 0;
x1 = pi/2;

fI = integral(f, x0, x1)


%%
%Regla del trapecio
clear all
clc

f = @(x) 6 + 3*cos(x);
x0 = 0;
x1 = pi/2;

h = x1 - x0;

fI = h/2*(f(x0) + f(x1))

ErrorRelativo(fI)

%%
%Regla del trapecio multiple con n = 2 y n = 4
clear all
clc

f = @(x) 6 + 3*cos(x);
x0 = 0;
x1 = pi/2;

fI2 = TrapecioCompuesto(f, x0, x1, 2);

fI4 = TrapecioCompuesto(f, x0, x1, 4);

ErrorRelativo(fI2)
ErrorRelativo(fI4)


%%
%Regla de simpson 1/3
clear all
clc

f = @(x) 6 + 3*cos(x);
x0 = 0;
x2 = pi/2;
h = (x2 - x0)/2
x1 = x0 + h;
fI = h/3*(f(x0) + 4*f(x1) + f(x2))
ErrorRelativo(fI)


%%
%Regla de simpson multiple 1/3 con n = 4
clear all
clc

f = @(x) 6 + 3*cos(x);
x0 = 0;
x2 = pi/2;


fi = SimpsonCompuesto(f, x0, x2, 4, 1/3)

ErrorRelativo(fi)

%%
%Regla de simpson 3/8
clear all
clc

f = @(x) 6 + 3*cos(x);
x0 = 0;
x2 = pi/2;
h = (x2 - x0)/2;
x1 = x0 + h;
fI = h*3/8*(f(x0) + 4*f(x1) + f(x2))

ErrorRelativo(fI)
%%
%Regla de simpson multiple 3/8 con n = 6
clear all
clc

f = @(x) 6 + 3*cos(x);
x0 = 0;
x2 = pi/2;


fI = SimpsonCompuesto(f, x0, x2, 6, 3/8)

ErrorRelativo(fI)

%% FUCIONES UTILIZADAS
function fI = TrapecioCompuesto(f, a, b, n)
    h = (b - a) / n;
    fI = f(a) + f(b);
    for i = 1:n-1
        X = a + i * h;
        fI = fI + 2 * f(X);
    end
    fI = h * fI / 2;
end

function fI = SimpsonCompuesto(f, a, b, n, fac)
    h = (b - a) / n;
    XI0 = f(a) + f(b);
    XI1 = 0; % Suma de f(x_2i-1)
    XI2 = 0; % Suma de f(x_2i)
    for i = 1:n-1
        X = a + i * h;
        if mod(i, 2) == 0
            XI2 = XI2 + f(X);
        else
            XI1 = XI1 + f(X);
        end
    end
    fI = h * (XI0 + 2 * XI2 + 4 * XI1) * fac;
end

function E = ErrorRelativo(aprox)
    E = ((abs(12.42478 - aprox) )/12.42478)*100;
    return
end