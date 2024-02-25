clc
syms X
sympref('FloatingPointOutput',true)

n = 5;
x = [1 ; 2 ; 3 ; 5 ; 7 ; 8];
y = [3 ; 6 ; 19 ; 99 ; 291 ; 444];


[a,b,c,d] = SplineCubico(n,x,y);
a
b
c
d

figure; hold on
for i = 1.0 :1: n 
        X = x(i):0.001:x(i+1);
        si = a(i, 1) + b(i, 1)*(X-x(i, 1)) + c(i, 1)*(X-x(i, 1)).^2 + d(i, 1)*(X-x(i, 1)).^3
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(x(i))+","+num2str(x(i+1))+"]"));
        xlim([0 10])
        %ylim([-0.4 0.8])
        legend();
end


function [a,b,c,d] = SplineCubico(n, x, y)
    h = zeros(n,1);
    alpha = zeros(n,1);

    a = y(1:n);
    c = zeros(n+1,1);
    b = zeros(n,1);
    d = zeros(n,1);

    I = zeros(n,1);
    I(1,1) = 1;
    U = zeros(n,1);
    Z = zeros(n,1);
    

    for i = 1.0 :1: n 
        h(i, 1) = x(i+1,1) - x(i,1);
    end
    


    for i = 2.0 :1: n
        alpha(i, 1) = (3/h(i, 1))*(y(i+1,1) - y(i,1)) - (3/h(i-1, 1))*(y(i,1) - y(i-1,1));
    end

    for i = 2.0 :1: n
        I(i,1) = 2*(x(i+1,1) - x(i-1,1)) - h(i-1, 1)*U(i-1, 1);
        U(i,1) = (h(i, 1))/(I(i, 1));
        Z(i,1) = ( alpha(i,1) - ( h(i-1,1) * Z(i-1,1) ))/I(i, 1);
    end

    I(n+1,1) = 1;
    U(n+1,1) = 0;
    Z(n+1,1) = 0;

    for j = n : -1 : 1
        c(j,1) = Z(j,1) - U(j,1)*c(j+1,1);
        b(j,1) = ( ( y(j+1,1) - y(j,1) )/(h(j,1)) ) - ( h(j,1)*(c(j+1,1) + 2*c(j,1) )/3 );
        d(j,1) = ( c(j+1,1) - c(j,1) )/(3*h(j,1));
    end

    c = c(1:n);

    return 
end