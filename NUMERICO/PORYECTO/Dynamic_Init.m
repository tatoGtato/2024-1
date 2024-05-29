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
T = 2;
dt = 0.1;
rc = 0.5;
%Antes de empezar el algoritmo se tienen que calcular el tamaño de la caja
%donde estaran las particulas (boundary conditions) 

%Se empiza suponiendo una masa de 0.0002 (creo q esto esta mal xd pero mientras) 
m = 1;
v = m/d;

box_ = v^(1/3);
%Con el volumen se saca BOX

%%
clc

molecularDynamics(T, dt, n, temp, box_, rc)

%% FUNCIONES UTILIZADAS

%ESTA FUNCION ES COMO EL MAIN DEL PROGRAMA, ACA SE EJECUTARAN TODAS
%LAS SUBRUTINAS
% function molecularDynamics(T, dt, n, temp, box, rc)
%    rc = min(rc, box);
% 
%    [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = dynamic_init(n, temp, dt, box);
%    subplot(1, 2, 1)
%    scatter3(Xx, Xy, Xz)
%    t = 0;
% 
%    while (t < T)
%        [Fx, Fy, Fz, en] = force(n, box, rc, Xx, Xy, Xz);
%        [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = integrate(n, dt, Fx, Fy, Fz, Vx, Vy, Vz, Xx, Xy, Xz, xmx, xmy, xmz);
%        t = t + dt
%    end 
% 
%    subplot(1, 2, 2)
%    scatter3(Xx, Xy, Xz)
%    Xx, Xy, Xz
% end
function molecularDynamics(T, dt, n, temp, box_, rc)
    rc = min(rc, box_);

    [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = dynamic_init(n, temp, dt, box_);
    subplot(1, 2, 1)
    scatter3(Xx, Xy, Xz)
    t = 0;

    while (t < T)
        [Fx, Fy, Fz, en] = force(n, box_, rc, Xx, Xy, Xz);
        [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = integrate(n, dt, Fx, Fy, Fz, Vx, Vy, Vz, Xx, Xy, Xz, xmx, xmy, xmz);
        t = t + dt

        % Llamada a la función SAMPLE para realizar el muestreo
        Switch = 1; % Establece el switch en 1 para indicar que es una llamada de muestreo
        Is = round(t / dt); % Calcula el número total de pasos de tiempo desde el inicio de la simulación
        En = en; % Energía total (potencial + cinética)
        Vir = sum(Fx.*Xx + Fy.*Xy + Fz.*Xz); % Virial total
        Enk = sum(0.5 * (Vx.^2 + Vy.^2 + Vz.^2)); % Energía cinética total
        Delt = dt; % Paso de tiempo de la simulación MD
        sample_2(Switch, Is, En, Vir, Enk, Delt); % Llama a la función SAMPLE
    end 

    subplot(1, 2, 2)
    scatter3(Xx, Xy, Xz)
    disp([Xx, Xy, Xz]); % Muestra las posiciones finales
end


%ESTA FUNCION NOS DA LA VELOCIDAD Y LAS POSICIONES DE LAS PARTICULAS
function [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = dynamic_init(n, temp, dt, box_)
    sumVx = 0;
    sumVy = 0;
    sumVz = 0;

    sumVx2 = 0;
    sumVy2 = 0;
    sumVz2 = 0;

    % Xx = zeros(n);  
    % Xy = zeros(n);  
    % Xz = zeros(n);

    Vx = zeros(n); 
    Vy = zeros(n);
    Vz = zeros(n);

    xmx = zeros(n);
    xmy = zeros(n);
    xmz = zeros(n);
    
    [Xx, Xy, Xz] = lattice_pos(n, box_);

    for i = 1 : n       %Seponen las particulas en una lattice con velocidad aleatoria 
        Vx(i) = randi(100) - 0.5;
        Vy(i) = randi(100) - 0.5;
        Vz(i) = randi(100) - 0.5;

        sumVx = sumVx + Vx(i);
        sumVy = sumVy + Vy(i);
        sumVz = sumVz + Vz(i);

        sumVx2 = sumVx2 + Vx(i)^2;
        sumVy2 = sumVy2 + Vy(i)^2;
        sumVz2 = sumVz2 + Vz(i)^2;
    end 

    sumVx = sumVx/n;
    sumVy = sumVy/n;
    sumVz = sumVz/n;

    sumVx2 = sumVx2/n;
    sumVy2 = sumVy2/n;
    sumVz2 = sumVz2/n;

    fsx = sqrt(3*temp/sumVx2);
    fsy = sqrt(3*temp/sumVy2);
    fsz = sqrt(3*temp/sumVz2);

    for i = 1 : n
        Vx(i) = (Vx(i) - sumVx)*fsx;
        Vy(i) = (Vy(i) - sumVy)*fsy;
        Vz(i) = (Vz(i) - sumVz)*fsz;

        xmx(i) = Xx(i) - Vx(i)*dt;
        xmy(i) = Xy(i) - Vy(i)*dt;
        xmz(i) = Xz(i) - Vz(i)*dt;
    end 
    return
end 

%ESTA ACTUA COMO FUNCION HELPER DE DYNAMIC_INIT() Y NOS DA UNA POSICION EN 
function [X, Y, Z] = lattice_pos(npart, box_)
    %LATTICE Place `npart` particles on a simple cubic lattice with density `rho`
    %
    % [X, Y, Z] = lattice(npart, box)
    % npart - number of particles
    % box - size of the box

    % Calculate the number of particles per side of the cubic lattice
    n = floor(npart^(1/3)) + 1;
    if n == 0
        n = 1;
    end

    % Calculate the distance between particles
    del = box_ / double(n);

    % Initialize particle coordinates
    X = zeros(1, npart);
    Y = zeros(1, npart);
    Z = zeros(1, npart);

    itel = 0;
    dx = -del;
    
    for i = 1:n
        dx = dx + del;
        dy = -del;
        for j = 1:n
            dy = dy + del;
            dz = -del;
            for k = 1:n
                dz = dz + del;
                if itel < npart
                    itel = itel + 1;
                    X(itel) = dx;
                    Y(itel) = dy;
                    Z(itel) = dz;
                end
            end
        end
    end
    
    fprintf('Initialisation on lattice: \n%d particles placed on a lattice\n', itel);
end

%ESTA FUNCION NOS DA A FUERZA Y LA ENERGIA DEL SISTEMA
function [Fx, Fy, Fz, en] = force(n, box_, rc, Xx, Xy, Xz)
    RC2 = rc*rc;

    en = 0;   

    Fx = zeros(n);
    Fy = zeros(n);
    Fz = zeros(n);


    for i = 1 : n-1
        for j = i+1 : n
            dx = Xx(i) - Xx(j);
            dy = Xy(i) - Xy(j);
            dz = Xz(i) - Xz(j);

            dx = dx-box_*round(dx/box_);
            dy = dy-box_*round(dy/box_);
            dz = dz-box_*round(dz/box_);

            r2 = dx^2 + dy^2 + dz^2;
 
            if (r2 < RC2)  
                r2i = 1/r2;
                r6i = r2i^3;
                ff = 48 * r2i * r6i * (r6i - 0.5);

                Fx(i) = Fx(i) + ff*dx;
                Fy(i) = Fy(i) + ff*dy;
                Fz(i) = Fz(i) + ff*dz;

                Fx(j) = Fx(j) - ff*dx;
                Fy(j) = Fy(j) - ff*dy;
                Fz(j) = Fz(j) - ff*dz;

                en = en + 4 * r6i * (r6i -1) - 4*((1/rc^12) - (1/rc^6));
            end
        end
    end
end

% %ESTA FUNCION INTEGRA LAS ECUACIONES DE MOVIMIENTO
function [XxNew, XyNew, XzNew, VxNew, VyNew, VzNew, xmx, xmy, xmz] = integrate(n, dt, Fx, Fy, Fz, Vx, Vy, Vz, Xx, Xy, Xz, xmx, xmy, xmz)
    v2 = 0;

    VxNew = zeros(n);
    VyNew = zeros(n);
    VzNew = zeros(n);

    XxNew = zeros(n);
    XyNew = zeros(n);
    XzNew = zeros(n);

    for i = 1:n
        vxt = Vx(i);
        vyt = Vy(i);
        vzt = Vz(i);

        VxNew(i) = Vx(i) + dt*Fx(i);
        VyNew(i) = Vy(i) + dt*Fy(i);
        VzNew(i) = Vz(i) + dt*Fz(i);

        XxNew(i) = Xx(i) + dt*VxNew(i);
        XyNew(i) = Xy(i) + dt*VyNew(i);
        XzNew(i) = Xz(i) + dt*VzNew(i);

        xmx(i) = XxNew(i);
        xmy(i) = XyNew(i);
        xmz(i) = XzNew(i);

        v2 = v2 + (Vx(i) + vxt)^(2/4) + (Vy(i) + vyt)^(2/4) + (Vz(i) + vzt)^(2/4);
    end
end 
