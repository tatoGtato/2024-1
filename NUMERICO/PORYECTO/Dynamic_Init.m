%% DEFINICION DE VARIABLES
% d = densidad
% temp = Temperatura
% n = # de particulas
% T = Tiempo maximo
% dt = Paso del tiempo
% m = Masa
% v = volumen
% box = Longitud de uno de los vertices de la caja
% x = Vector de las posiciones de las particuals
% v = Vector de las velocidades de las particulas
% f = Vector de fuerzas del sistema
% en = Energia del sistema
% etot = Energia total por particula
% 

%%
clear all 
clc 

d = 0.8442;
temp = 0.728;
n = 100;
T = 100;
dt = 0.01;

%Antes de empezar el algoritmo se tienen que calcular el tamaño de la caja
%donde estaran las particulas (boundary conditions) 

%Se empiza suponiendo una masa de 0.0002 (creo q esto esta mal xd pero mientras) 
m = 0.0002;
v = m/d;

box = v^(1/3);
%Con el volumen se saca BOX

%%

molecularDynamics(T, dt, n, d, temp, box)

%% FUNCIONES UTILIZADAS

%ESTA FUNCION ES COMO EL MAIN DEL PROGRAMA, ACA SE EJECUTARAN TODAS
%LAS SUBRUTINAS
function molecularDynamics(T, dt, n, d, temp, box)

   [x,xm, v] = dynamic_init(n, temp, dt);
   t = 0;


   while (t < T)
       [f, en] = force(n, x, box);
       integrate(f, en, n, x, xm, dt);
       t = t + dt;
   end 
end 


%ESTA FUNCION NOS DA LA VELOCIDAD Y LAS POSICIONES DE LAS PARTICULAS
function [x,xm] = dynamic_init(n, temp, dt)
    sumv = 0;  %Se inicializa la suma de velocidades
    sumv2 = 0;
    x = zeros(n);  %Vector de posiciones
    v = zeros(n);   %Vector de velocidades
    xm = zeros(n);
    
    for i = 1 : n       %Seponen las particulas en una lattice con velocidad aleatoria
        x(i) = lattice_pos(i, n);
        v(i) = randi(100) - 0.5;
        sumv = sumv + v(i);
        sumv2 = sumv2 + v(i)^2;
    end 

    sumv = sumv/n;
    sumv2 = sumv2/n;

    fs = sqrt(3*temp/sumv2);

    for i = 1 : n
        v(i) = (v(i) - sumv)*fs;
        xm(i) = x(i) - v(i)*dt;
    end 
    return
end 

%ESTA ACTUA COMO FUNCION HELPER DE DYNAMIC_INIT() Y NOS DA UNA POSICION EN 
%UN LATTICE 
function [x, y] = lattice_pos(i, n)
    % Convert the linear index to 2D coordinates
    y = ceil(i / n);       % Row index
    x = mod(i-1, n) + 1;   % Column index
end


%ESTA FUNCION NOS DA A FUERZA Y LA ENERGIA DEL SISTEMA
function [f, en] = force(n, x, box)
    en = 0;   
    f = zeros(n);

    for i = 1 : n-1
        for j = i+1 : n
            xr = x(i) - x(j);
            xr = xr-box*nint(xr/box);
            r2 = xr^2;

            if (r2 < rc*2)
                r2i = 1/r2;
                r6i = r2i^3;
                ff = 48 * r2i * r6i * (r6i - 0.5);
                f(i) = f(i) + ff*xr;
                f(j) = f(j) - ff*xr;
                en = en + 4 * r6i * (r6i -1) - 4*((1/rc^12) - (1/rc^6));
            end
        end
    end
end

%ESTA FUNCION INTEGRA LAS ECUACIONES DE MOVIMIENTO
function [etot, tempI] = integrate(f, en, n, x, xm, dt)
    sumv = 0;
    sumv2 = 0;

    for i = 1 :n 
        xx = 2*x(i) - xm(i) + dt^2 * f(i);
        vi = (xx-xm(i))/(2*dt);
        sumv = sumv + vi;
        sumv2 =sumv2 + vi^2;
        xm(i) = x(i);
        x(i) = xx;
    end 
    tempI = sumv2/(3*n);
    etot = (en + 0.5*sumv2)/n;
end 





