%%
%PUNTO 1 - GRAFICAR LA NUBE DE PUNTOS

clear all
syms X
clc

xArriba = [16;17;18;19;20;20.5;21;22;23;24;25;26;27;27.6;27.8;29;30;31;32];
yArriba = [11.5;12.3;12.7;13.1;13.4;13.5;14.6;15.5;15.9;16.1;16;15.7;15;14;13.4;13.1;12.8;12.3;11.5];

xAbajo = [16;17;18;19;20;20.5;21;22;23;24;25;26;27;27.5;28;29;30;31;32];
yAbajpo = [11.5;10.8;10.3;9.9;9.6;9.5;8.8;8.4;8;7.9;8;8.3;8.7;9.7;9.8;10;10.5;10.9;11.5];

xtotal = [];
ytotal = [];

hold on 
plot(xAbajo,yAbajpo, 'o')
plot(xArriba,yArriba, 'o')
xlim([6 40])
ylim([6 18])

%%
% Punto 2 - Graficar el contorno generado por los splines

%ARRIBA
%ALA 1
xArribaAla1 = [16;17;18;19;20;20.5];
yArribaAla1 = [11.5;12.3;12.7;13.1;13.4;13.5];

%CIRCULO
xArribaCirc = [20.5;21;22;23;24;25;26;27;27.6;27.8];
yArribaCirc = [13.5;14.6;15.5;15.9;16.1;16;15.7;15;14;13.4];

%ALA 2
xArribaAla2 = [27.8;29;30;31;32];
yArribaAla2 = [13.4;13.1;12.8;12.3;11.5];


%ABAJO
%ALA 1
xAbajoAla1 = [16;17;18;19;20;20.5];
yAbajoAla1 = [11.5;10.8;10.3;9.9;9.6;9.5];


%CIRCULO
xAbajoCirc = [20.5;21;22;23;24;25;26;27;27.5];
yAbajoCirc = [9.5;8.8;8.4;8;7.9;8;8.3;8.7;9.7];


%ALA 2
xAbajoAla2 = [27.5;28;29;30;31;32];
yAbajoAla2 = [9.7;9.8;10;10.5;10.9;11.5];


%SPLINES PARA ARRIBA
%ala1
nArala1 = size(xArribaAla1,1) - 1;
[aArala1, bArala1, cArala1, dArala1] = SplineCubico(nArala1, xArribaAla1, yArribaAla1);
%circulo
nArCirc = size(xArribaCirc,1) - 1;
[aArCirc, bArCirc, cArCirc, dArCirc] = SplineCubico(nArCirc, xArribaCirc, yArribaCirc);
%ala2
nArala2 = size(xArribaAla2,1) - 1;
[aArala2, bArala2, cArala2, dArala2] = SplineCubico(nArala2, xArribaAla2, yArribaAla2);

%SPLINE PARA ABAJO
%ala1
nAbala1 = size(xAbajoAla1,1) - 1;
[aAbala1, bAbala1, cAbala1, dAbala1] = SplineCubico(nAbala1, xAbajoAla1, yAbajoAla1);
%circulo
nAbCirc = size(xAbajoCirc,1) - 1;
[aAbCirc, bAbCirc, cAbCirc, dAbCirc] = SplineCubico(nAbCirc, xAbajoCirc, yAbajoCirc);
%ala2
nAbala2 = size(xAbajoAla2,1) - 1;
[aAbala2, bAbala2, cAbala2, dAbala2] = SplineCubico(nAbala2, xAbajoAla2, yAbajoAla2)


%Para plotear 
figure; hold on
%Arriba
%ala1
for i = 1.0 :1: nArala1
        X = xArribaAla1(i):0.001:xArribaAla1(i+1);  
        si = aArala1(i, 1) + bArala1(i, 1)*(X-xArribaAla1(i, 1)) + cArala1(i, 1)*(X-xArribaAla1(i, 1)).^2 + dArala1(i, 1)*(X-xArribaAla1(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(xArribaAla1(i))+","+num2str(xArribaAla1(i+1))+"]"));
        xlim([6 40])
        ylim([6 18])
        legend();
end
% %circ
for i = 1.0 :1: nArCirc
        X = xArribaCirc(i):0.001:xArribaCirc(i+1);  
        si = aArCirc(i, 1) + bArCirc(i, 1)*(X-xArribaCirc(i, 1)) + cArCirc(i, 1)*(X-xArribaCirc(i, 1)).^2 + dArCirc(i, 1)*(X-xArribaCirc(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(xArribaCirc(i))+","+num2str(xArribaCirc(i+1))+"]"));
        xlim([6 40])
        ylim([6 18])
        legend();
end
% %ala2
for i = 1.0 :1: nArala2
        X = xArribaAla2(i):0.001:xArribaAla2(i+1);  
        si = aArala2(i, 1) + bArala2(i, 1)*(X-xArribaAla2(i, 1)) + cArala2(i, 1)*(X-xArribaAla2(i, 1)).^2 + dArala2(i, 1)*(X-xArribaAla2(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(xArribaAla2(i))+","+num2str(xArribaAla2(i+1))+"]"));
        xlim([6 40])
        ylim([6 18])
        legend();
end
% %Abajo
% %ala1
for i = 1.0 :1: nAbala1
        X = xAbajoAla1(i):0.001:xAbajoAla1(i+1);  
        si = aAbala1(i, 1) + bAbala1(i, 1)*(X-xAbajoAla1(i, 1)) + cAbala1(i, 1)*(X-xAbajoAla1(i, 1)).^2 + dAbala1(i, 1)*(X-xAbajoAla1(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(xAbajoAla1(i))+","+num2str(xAbajoAla1(i+1))+"]"));
        xlim([6 40])
        ylim([6 18])
        legend();
end
%circulo
for i = 1.0 :1: nAbCirc
        X = xAbajoCirc(i):0.001:xAbajoCirc(i+1);  
        si = aAbCirc(i, 1) + bAbCirc(i, 1)*(X-xAbajoCirc(i, 1)) + cAbCirc(i, 1)*(X-xAbajoCirc(i, 1)).^2 + dAbCirc(i, 1)*(X-xAbajoCirc(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(xAbajoCirc(i))+","+num2str(xAbajoCirc(i+1))+"]"));
        xlim([6 40])
        ylim([6 18])
        legend();
end
%ala2
for i = 1.0 :1: nAbala2
        X = xAbajoAla2(i):0.001:xAbajoAla2(i+1);  
        si = aAbala2(i, 1) + bAbala2(i, 1)*(X-xAbajoAla2(i, 1)) + cAbala2(i, 1)*(X-xAbajoAla2(i, 1)).^2 + dAbala2(i, 1)*(X-xAbajoAla2(i, 1)).^3;
        plot(X,si,'DisplayName',strcat('S=',"["+num2str(xAbajoAla2(i))+","+num2str(xAbajoAla2(i+1))+"]"));
        xlim([6 40])
        ylim([6 18])
        legend();
end


%%
%Punto 3 - Utilizar simpson para aproximar el area  

%Area de la curva de arriba
areaArribaala1 = 0;
areaArribacirc = 0;
areaArribaala2 = 0;


for i = 1.0 :1: nArala1
    si = @(X) aArala1(i, 1) + bArala1(i, 1)*(X-xArribaAla1(i, 1)) + cArala1(i, 1)*(X-xArribaAla1(i, 1)).^2 + dArala1(i, 1)*(X-xArribaAla1(i, 1)).^3;
    areaArribaala1 = areaArribaala1 + SimpsonCompuesto(si, xArribaAla1(i), xArribaAla1(i+1), 4);
end
for i = 1.0 :1: nArCirc
    si = @(X) aArCirc(i, 1) + bArCirc(i, 1)*(X-xArribaCirc(i, 1)) + cArCirc(i, 1)*(X-xArribaCirc(i, 1)).^2 + dArCirc(i, 1)*(X-xArribaCirc(i, 1)).^3;
    areaArribacirc = areaArribacirc + SimpsonCompuesto(si, xArribaCirc(i), xArribaCirc(i+1), 4);
end
for i = 1.0 :1: nArala2
    si = @(X) aArala2(i, 1) + bArala2(i, 1)*(X-xArribaAla2(i, 1)) + cArala2(i, 1)*(X-xArribaAla2(i, 1)).^2 + dArala2(i, 1)*(X-xArribaAla2(i, 1)).^3;
    areaArribaala2 = areaArribaala2 + SimpsonCompuesto(si, xArribaAla2(i), xArribaAla2(i+1), 4);
end

areaArribaTotal = areaArribaala1 + areaArribacirc + areaArribaala2;

%Area de la curva de abajo
areaAbajoala1 = 0;
areaAbajocirc = 0;
areaAbajoala2 = 0;

for i = 1.0 :1: nAbala1
    si = @(X) aAbala1(i, 1) + bAbala1(i, 1)*(X-xAbajoAla1(i, 1)) + cAbala1(i, 1)*(X-xAbajoAla1(i, 1)).^2 + dAbala1(i, 1)*(X-xAbajoAla1(i, 1)).^3;
    areaAbajoala1 = areaAbajoala1 + SimpsonCompuesto(si, xAbajoAla1(i), xAbajoAla1(i+1), 4);
end
for i = 1.0 :1: nAbCirc
    si = @(X) aAbCirc(i, 1) + bAbCirc(i, 1)*(X-xAbajoCirc(i, 1)) + cAbCirc(i, 1)*(X-xAbajoCirc(i, 1)).^2 + dAbCirc(i, 1)*(X-xAbajoCirc(i, 1)).^3;
    areaAbajocirc = areaAbajocirc + SimpsonCompuesto(si, xAbajoCirc(i), xAbajoCirc(i+1), 4);   
end
for i = 1.0 :1: nAbala2
    si = @(X) aAbala2(i, 1) + bAbala2(i, 1)*(X-xAbajoAla2(i, 1)) + cAbala2(i, 1)*(X-xAbajoAla2(i, 1)).^2 + dAbala2(i, 1)*(X-xAbajoAla2(i, 1)).^3;
    areaAbajoala2 = areaAbajoala2 + SimpsonCompuesto(si, xAbajoAla2(i), xAbajoAla2(i+1), 4);
end

areaAbajoTotal = areaAbajoala1 + areaAbajocirc + areaAbajoala2

areaOvni = areaArribaTotal - areaAbajoTotal
areaImagen = 26*46

porcentaje = (areaOvni/areaImagen)*100

%%
%FUNCIONES USADAS
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

function fI = SimpsonCompuesto(f, a, b, n)
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
    fI = h * (XI0 + 2 * XI2 + 4 * XI1) * 1/3;
end


