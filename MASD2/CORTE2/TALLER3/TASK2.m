%%
%Aproximada
clc
clear all 

K = 1;
dx = pi/20;
dt = (pi^2/800);
muuuu = (K*dt)/(dx^2);
T = (3*pi^2/80);
Uo = zeros(1,ceil(T/dt)+1);


%Construimos Uo
i = 1;
for x = 0.0: dx :pi
    if (x < pi/2)
        Uo(i) = x;

    elseif (x == pi/2)
        Uo(i) = pi/2;

    else
        Uo(i) = pi-x;
    end
    i = i +1;
end

%Construimos la aproximacion

for j = 0: dt : T
    U = Uo.*(1 - 4*(muuuu)*sin((K*dx)/2)^2);
    Uo = U;
end



figure
X = 0: dt: T;
Y = Uo;
size(X)
size(Y)
stem(X,Y)


%%
%EXACTA
syms x t 

bk = @(k) (4*(-1)^((k+1)/2))/pi*k^2;

Dx = pi/20;

Dt = pi/40;

U = @(k) bk(k)*sin(k*x)*exp(-k^2*t)

Uf = U(1)

for k = 2: 1 : 100
    Uf = Uf + U(k);
end

figure
Y = U;
X = linspace(0,1,21);
stem(X,Y)
