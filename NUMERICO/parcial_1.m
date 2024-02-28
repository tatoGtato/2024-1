%%
%PUNTO 1

clc
syms X
sympref('FloatingPointOutput',true)

x = [1.74 ; 2.72 ; 3.72; 4.09 ; 4.32 ; 4.70 ; 5.00 ; 6.00 ; 6.53 ; 6.70 ; 7.29; 8.06; 10.02; 11.12 ];
y = [-5.3 ; -10.08 ;-21.8 ; -32.0 ; -35.8 ; -36.7 ; -36.7 ; -33.2 ; -15.7 ; -10.0 ; 13.7 ; 32.2 ; 24.0 ; 6.9]
size(x)
size(y)

n = 13;

[a, b, c, d] = SplineCubico(n,x, y)

figure; hold on
for i = 1.0 :1: n
        X = x(i):0.001:x(i+1);  
        si = a(i, 1) + b(i, 1)*(X-x(i, 1)) + c(i, 1)*(X-x(i, 1)).^2 + d(i, 1)*(X-x(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(x(i))+","+num2str(x(i+1))+"]"));
        xlim([0 12])
        ylim([-40 40])
        legend();
end

%%
%PUNTO 2 (Proceso)
clc
syms X
sympref('FloatingPointOutput',true)

x = [1.74 ; 2.72 ; 3.72; 4.09 ; 4.32 ; 4.70 ; 5.00 ; 6.00 ; 6.53 ; 6.70 ; 7.29; 8.06; 10.02; 11.12 ];
y = [-5.3 ; -10.08 ;-21.8 ; -32.0 ; -35.8 ; -36.7 ; -36.7 ; -33.2 ; -15.7 ; -10.0 ; 13.7 ; 32.2 ; 24.0 ; 6.9]
size(x)
size(y)

n = 13;

[a, b, c, d] = SplineCubico(n,x, y)


figure; hold on
for i = 1.0 :1: n
        if (x(i) == 6.7 && x(i+1) == 7.29 )
            si = a(i, 1) + b(i, 1)*(X-x(i, 1)) + c(i, 1)*(X-x(i, 1)).^2 + d(i, 1)*(X-x(i, 1)).^3
        end
end

f = @(z) 35.2137*z + 22.1808*(z - 6.7000)^2 - 23.3580*(z - 6.7000)^3 - 245.9318;

tol = 0.0005;
No = 20;

p0 = 7.28;
p1 = 6.8;
%CON FALSA POSICION
[p,i] = falsaPosicion(f, p0, p1, tol, No)

%%
%PUNTO 3 (PLOTEAR LA GRAFICA CON EL PUNTO EN 0)

clc
syms X
sympref('FloatingPointOutput',true)

x = [1.74 ; 2.72 ; 3.72; 4.09 ; 4.32 ; 4.70 ; 5.00 ; 6.00 ; 6.53 ; 6.70 ; 7.29; 8.06; 10.02; 11.12 ];
y = [-5.3 ; -10.08 ;-21.8 ; -32.0 ; -35.8 ; -36.7 ; -36.7 ; -33.2 ; -15.7 ; -10.0 ; 13.7 ; 32.2 ; 24.0 ; 6.9]
size(x)
size(y)

n = 13;

[a, b, c, d] = SplineCubico(n,x, y)

f = @(z) 35.2137*z + 22.1808*(z - 6.7000)^2 - 23.3580*(z - 6.7000)^3 - 245.9318;

punto0x = 6.9542;
punto0y = f(6.9542);

figure; hold on
for i = 1.0 :1: n
        X = x(i):0.001:x(i+1);  
        si = a(i, 1) + b(i, 1)*(X-x(i, 1)) + c(i, 1)*(X-x(i, 1)).^2 + d(i, 1)*(X-x(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(x(i))+","+num2str(x(i+1))+"]"));
        xlim([0 12])
        ylim([-40 40])
        legend();
end
plot(punto0x,punto0y,'o')


%%
%SPLINE
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
    

    for i = 1.0 :1: n tol = 0.0005;
No = 20;
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

%%
%FALSA-POSICION

function [p,i] = falsaPosicion(f, p0, p1, tol, No)

    i = 2;
    q0 = f(p0);
    q1 = f(p1);
    
    while i <= No
        p = p1 - ( (q1*(p1-p0))/(q1-q0) );
        if abs(p - p1) < tol 
            return
        end
        i = i+1;
        q = f(p);
        if (q * q1 < 0)
            p0 = p1;
            q0 = q1;
        end
        p1 = p;
        q1 = q;
        
    end
end



