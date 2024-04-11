clc
clear all

%EJ 1
n = 4;
a = [4 1 -1 1 ;...
     1 4 -1 -1; ...
     -1 -1 5 1; ...
     1 -1 1 3];
b = [-2 ; -1 ; 0 ; 1];
Xo = [0 ; 0 ; 0 ; 0];
TOL = 10*10^(-6);
nMax = 30;

[x, k] = jacop(n, a, b, Xo,TOL, nMax)

%% EJ 2
clc
clear all

%EJ 1
n = 3;
a = [1 2 -2  ;...
     1 1 1; ...
     2 2 1;];
b = [7 ; 2 ; 5];
Xo = [0 ; 0 ; 0];
TOL = 10*10^(-6);
nMax = 30;

[x, k] = jacop(n, a, b, Xo,TOL, nMax)
D = diag(diag(a))
L = tril(a, -1)
U = triu(a, 1)

Tj = -inv(D).*(L+U)

eigj = eig(Tj)

%%

function [x, k] = jacop(n, a, b, Xo, TOL, nMax)
    k = 1;
    x = zeros(n,1);
    while (k <= nMax)
        for i = 1 : n
            x(i) = 1/a(i,i) * (b(i) - suma(i,n,a,Xo)); %HACER SUMA BIEN
        end
        if (norm(x - Xo)/norm(x) < TOL)
            return 
        end
        k = k + 1;
        for i = 1 : n
            Xo(i) = x(i);
        end 
    end 
    "SE EXEDIO EL NUM DE ITERACIONES"
end 

function sum = suma(i, n, a, Xo)
    sum = 0;
    for j = 1 : n
        if (j ~= i)
            sum = sum + a(i,j)*Xo(j);
        end 
    end 
end 