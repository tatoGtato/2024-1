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

[x, k] = gausSeidel(n, a, b, Xo,TOL, nMax)

%%
function [x, k] = gausSeidel(n, a, b, Xo, TOL, nMax)
    k = 1;
    x = zeros(n,1);
    while (k <= nMax)
        for i = 1 : n
            x(i) = 1/a(i,i) * (b(i) - suma1(i,a,x) - suma2(i, n, a, Xo)); %HACER SUMA BIEN
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

function sum = suma1(i, a, x)
    sum = 0;
    for j = 1 : i-1   
        sum = sum + a(i,j)*x(j);
    end 
end 

function sum = suma2(i, n, a, Xo)
    sum = 0;
    for j = i+1 : n
        
            sum = sum + a(i,j)*Xo(j);
        
    end 
end 


