clc
clear all
syms X
sympref('FloatingPointOutput',true)

x=[0;20;80;100;130;150;170;200;235;250;270;300;320;350;380;400;435;450;460;490];
yA=[0;0;0.1;0.25;1;1.35;1.55;1.85;2;2.2;2.5; 2.8; 3  ; 3.4; 4  ; 4.3; 5;   5.5; 6;   6.5];
yN=[0, 0.1,0.35,0.4,0.5,0.6, 0.85,0.95,0.85,0.6, 0.45,0.35,0.2,0.1, 0.05,0  , 0,  -0.2,0  ,0];

y = zeros(20,1);

distancia = zeros(20,1);
aceleracion = zeros(20,1);

for i = 1: size(x,1)
    y(i) = sqrt( (yA(i))^2 + (yN(i))^2 );
end

plot(x,y, 'o')
xlim([-1 500])
ylim([-1 7])

n = 19;

[a, b, c, d] = SplineCubico(n,x, y);

%Para plotear 
% figure; hold on
% for i = 1.0 :1: n
%         X = x(i):0.001:x(i+1);  
%         si = a(i, 1) + b(i, 1)*(X-x(i, 1)) + c(i, 1)*(X-x(i, 1)).^2 + d(i, 1)*(X-x(i, 1)).^3
%         plot(X,si,'DisplayName',strcat('S=',"["+num2str(x(i))+","+num2str(x(i+1))+"]"));
%         xlim([0 500])
%         ylim([-1 7])
%         legend();
% end

%Para simpson
for i = 1.0 :1: n 
     si = @(X) a(i, 1) + b(i, 1)*(X-x(i, 1)) + c(i, 1)*(X-x(i, 1)).^2 + d(i, 1)*(X-x(i, 1)).^3;
     distancia(i) = SimpsonCompuesto(si, x(i), x(i+1), n, 1/3);
end

%Para derivacion 
for i = 1.0 :1: n 
     si = @(X) a(i, 1) + b(i, 1)*(X-x(i, 1)) + c(i, 1)*(X-x(i, 1)).^2 + d(i, 1)*(X-x(i, 1)).^3;
     aceleracion(i) = richardson_extrapolation_3puntos(si, x(i), 0.2);
end

%%



%%

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

function [derivada] = richardson_extrapolation_3puntos(f, x, h)
    %3 puntos
    derivada = (1/(2*h)) * (-3*f(x) + 4*f(x + h) - f(x + 2*h));

    % Extrapolación de Richardson
    f1 = richardson_aux(f, x, h, 1);
    f2 = richardson_aux(f, x, h, 2);

    % Corrección de la derivada utilizando extrapolación de Richardson
    derivada = (4*f1 - f2) / 3;
end

function [f] = richardson_aux(f, x, h, n)
    % Función auxiliar para la extrapolación de Richardson.
    f = (1/(2^n*h)) * (-f(x + h) + f(x - h));
end
