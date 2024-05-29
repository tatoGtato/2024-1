%% DEFINICION DE VARIABLES
clear all 
close all
clc 

d = 0.8442;
temp = 0.728;
n = 100;
T = 0.9;
dt = 0.01;
rc = 0.5;


% Calculo del tamaño de la caja
m = 1;
v = m/d;
rho = n/v
box_length = v^(1/3);

%% Inicialización de las variables de sample_2
Switch = 0; % Inicializa variables
Is = 0; % Se actualizará en el bucle principal
En = 0; % Se actualizará en el bucle principal
Vir = 0; % Se actualizará en el bucle principal
Enk = 0; % Se actualizará en el bucle principal
Delt = dt; % Paso de tiempo de la simulación MD

sample_2(Switch, Is, En, Vir, Enk, Delt, rho, temp, 0, 0, 0, 0, 0, 0);

%% Ejecución de la simulación
clc
molecularDynamics(T, dt, n, temp, box_length, rc, rho)

%% FUNCIONES UTILIZADAS

function molecularDynamics(T, dt, n, temp, box_length, rc, rho)
    rc = min(rc, box_length);

    [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = dynamic_init(n, temp, dt, box_length);
    subplot(3, 1, 1)
    scatter3(Xx, Xy, Xz)
    t = 0;
    ind = 0;

    while (t < T)
        [Fx, Fy, Fz, en] = force(n, box_length, rc, Xx, Xy, Xz);
        [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = integrate(n, dt, Fx, Fy, Fz, Vx, Vy, Vz, Xx, Xy, Xz, xmx, xmy, xmz, box_length);
        
        t = t + dt;
            
        % Cálculo de la energía cinética
        Enk = sum(0.5 * (Vx.^2 + Vy.^2 + Vz.^2));
    
        % Cálculo del virial
        Vir = sum(Fx.*Xx + Fy.*Xy + Fz.*Xz); 
    
        % Llamada a la función SAMPLE para realizar el muestreo
        Switch = 1; % Establece el switch en 1 para indicar que es una llamada de muestreo
        Is = round(t / dt); % Calcula el número total de pasos de tiempo desde el inicio de la simulación
        En = en + Enk; % Energía total (potencial + cinética)
        Delt = dt; % Paso de tiempo de la simulación MD
        sample_2(Switch, Is, En, Vir, Enk, Delt, rho, temp, Xx, Xy, Xz, Vx, Vy, Vz); % Llama a la función SAMPLE
        ind = ind + 1; 
    end 

    % Llamada a la función SAMPLE para escribir resultados al disco
    Switch = 2; % Establece el switch en 2 para indicar que se escriban los resultados
    r2t = sample_2(Switch, Is, En, Vir, Enk, Delt, rho, temp, Xx, Xy, Xz, Vx, Vy, Vz); % Llama a la función SAMPLE

    subplot(3, 1, 2)
    scatter3(Xx, Xy, Xz)

    subplot(3, 1, 3)
    plot(r2t)
end

function [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = dynamic_init(n, temp, dt, box_length)
    sumVx = 0;
    sumVy = 0;
    sumVz = 0;

    sumVx2 = 0;
    sumVy2 = 0;
    sumVz2 = 0;

    Vx = zeros(n, 1);
    Vy = zeros(n, 1);
    Vz = zeros(n, 1);

    xmx = zeros(n, 1);
    xmy = zeros(n, 1);
    xmz = zeros(n, 1);
    
    [Xx, Xy, Xz] = lattice_pos(n, box_length);

    for i = 1 : n       
        Vx(i) = rand - 0.5;
        Vy(i) = rand - 0.5;
        Vz(i) = rand - 0.5;

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
end

function [X, Y, Z] = lattice_pos(npart, box_length)
    n = floor(npart^(1/3)) + 1;
    if n == 0
        n = 1;
    end

    del = box_length / double(n);

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
                itel = itel + 1;
                if itel <= npart
                    X(itel) = dx;
                    Y(itel) = dy;
                    Z(itel) = dz;
                end
            end
        end
    end
end

function [Fx, Fy, Fz, en] = force(n, box_length, rc, Xx, Xy, Xz)
    Fx = zeros(n, 1);
    Fy = zeros(n, 1);
    Fz = zeros(n, 1);
    en = 0;
    
    rc2 = rc * rc;
    
    for i = 1:n-1
        for j = i+1:n
            dx = Xx(i) - Xx(j);
            dy = Xy(i) - Xy(j);
            dz = Xz(i) - Xz(j);
            
            % Consideración de condiciones periódicas de contorno
            if dx > 0.5 * box_length
                dx = dx - box_length;
            elseif dx < -0.5 * box_length
                dx = dx + box_length;
            end

            if dy > 0.5 * box_length
                dy = dy - box_length;
            elseif dy < -0.5 * box_length
                dy = dy + box_length;
            end

            if dz > 0.5 * box_length
                dz = dz - box_length;
            elseif dz < -0.5 * box_length
                dz = dz + box_length;
            end

            r2 = dx * dx + dy * dy + dz * dz;

            if r2 < rc2
                r2i = 1 / r2;
                r6i = r2i * r2i * r2i;
                ff = 48.0 * r2i * r6i * (r6i - 0.5);

                Fx(i) = Fx(i) + ff * dx;
                Fy(i) = Fy(i) + ff * dy;
                Fz(i) = Fz(i) + ff * dz;

                Fx(j) = Fx(j) - ff * dx;
                Fy(j) = Fy(j) - ff * dy;
                Fz(j) = Fz(j) - ff * dz;

                en = en + 4 * r6i * (r6i - 1.0);
            end
        end
    end
end

function [Xx, Xy, Xz, Vx, Vy, Vz, xmx, xmy, xmz] = integrate(n, dt, Fx, Fy, Fz, Vx, Vy, Vz, Xx, Xy, Xz, xmx, xmy, xmz, box_length)
    for i = 1 : n
        xx = 2 * Xx(i) - xmx(i) + Fx(i) * dt * dt;
        yy = 2 * Xy(i) - xmy(i) + Fy(i) * dt * dt;
        zz = 2 * Xz(i) - xmz(i) + Fz(i) * dt * dt;

        Vx(i) = (xx - xmx(i)) / (2 * dt);
        Vy(i) = (yy - xmy(i)) / (2 * dt);
        Vz(i) = (zz - xmz(i)) / (2 * dt);

        xmx(i) = Xx(i);
        xmy(i) = Xy(i);
        xmz(i) = Xz(i);

        Xx(i) = xx;
        Xy(i) = yy;
        Xz(i) = zz;

        if Xx(i) < 0
            Xx(i) = Xx(i) + box_length;
            xmx(i) = xmx(i) + box_length;
        elseif Xx(i) > box_length
            Xx(i) = Xx(i) - box_length;
            xmx(i) = xmx(i) - box_length;
        end

        if Xy(i) < 0
            Xy(i) = Xy(i) + box_length;
            xmy(i) = xmy(i) + box_length;
        elseif Xy(i) > box_length
            Xy(i) = Xy(i) - box_length;
            xmy(i) = xmy(i) - box_length;
        end

        if Xz(i) < 0
            Xz(i) = Xz(i) + box_length;
            xmz(i) = xmz(i) + box_length;
        elseif Xz(i) > box_length
            Xz(i) = Xz(i) - box_length;
            xmz(i) = xmz(i) - box_length;
        end
    end
end

function r2tRet = sample_2(Switch, Is, En, Vir, Enk, Delt, rho, temp, X, Y ,Z, VX, VY, VZ)
    % Sample averages:
    %   a) density, pressure, potential energy
    %   b) stress tensor correlation functions
    %   c) velocity autocorrelation function
    %   d) mean square displacement (conventional algorithm)
    %
    % Switch (input) = 1: sample averages
    %                = 0: initialize variables
    %                = 2: writes results to disk
    % Is      = total number of time steps since start simulation
    % En     (input) total energy (potential + kinetic)
    % Vir    (input) total virial
    % Delt   (input) time step md simulation

    persistent ngr tempav delg delgi g nhgr
    persistent vacf nvacf t0 vxt0 vyt0 vzt0 tvacf dtime r2t rx0 ry0 rz0 tt0 ttv0
    persistent sxyt sxzt syzt sxy0 sxz0 syz0 tstress0 dstresstime nstress sxy00 sxz00 syz00 tstress

    % Constants and Parameters
    NHIsmax = 250;
    TMAx = 1000;
    T0Max = 200;
    NPMax = 1000; % Assumption, define as needed
    PI = 3.141592653589793;
    HBOX = 10; % Assumption, define as needed
    BOX = 20; % Assumption, define as needed
    NPART = 100; % Assumption, define as needed
    RC2 = 9; % Assumption, define as needed
    ECUT = 0; % Assumption, define as needed
    NSAMP = 10; % Assumption, define as needed
    IGR = 100; % Assumption, define as needed
    NTVACF = 100; % Assumption, define as needed
    IT0 = 10; % Assumption, define as needed
    ITSTRESS0 = 10; % Assumption, define as needed
    TDIFMAX = 100; % Assumption, define as needed
    rc = BOX/2;
    
    if Switch == 1
        % Sample averages
        if NPART ~= 0
            enp = (En - Enk) / NPART;
            temp = 2 * Enk / (3 * NPART);
            vol = BOX^3;
            rho = NPART / vol; %densidad
            %press = rho * temp + Vir / (3 * vol);
            press = 16/3 * pi* rho^2 *((2/3)*((1/rc)^9) - (1/rc)^3); %presion 
        else
            rho = 0;
            enp = 0;
            press = 0; %presion
        end
        %fprintf(66, '%d %f %f %f\n', Is, temp, press, enp);
        %  SE CAMBIO EL ORDEN DE IGR Y IS
        if mod(IGR, Is) == 0
            % Sample radial distribution function and stress tensor
            ngr = ngr + 1;
            tempav = tempav + temp;
            sxy = 0;
            sxz = 0;
            syz = 0;
            for i = 1:NPART
                xi = X(i);
                yi = Y(i);
                zi = Z(i);
                for j = i + 1:NPART
                    dx = xi - X(j);
                    dy = yi - Y(j);
                    dz = zi - Z(j);
                    % Periodic boundary conditions
                    dx = dx - BOX * round(dx / BOX);
                    dy = dy - BOX * round(dy / BOX);
                    dz = dz - BOX * round(dz / BOX);
                    r2 = dx * dx + dy * dy + dz * dz;
                    if r2 <= (HBOX * HBOX)
                        r = sqrt(r2);
                        ig = floor(r * delgi);
                        g(ig) = g(ig) + 1;
                        % Stress tensor
                        if r2 <= RC2
                            r2i = 1 / r2;
                            r6i = r2i * r2i * r2i;
                            virij = 48 * (r6i * r6i - 0.5 * r6i);
                            fr = -virij * r2i;
                            fy = fr * dy;
                            fz = fr * dz;
                            sxy = sxy + dx * fy;
                            sxz = sxz + dx * fz;
                            syz = syz + dy * fz;
                        end
                    end
                end
                sxy = sxy + VX(i) * VY(i);
                sxz = sxz + VX(i) * VZ(i);
                syz = syz + VY(i) * VZ(i);
            end
            sxy00 = sxy00 + sxy * sxy;
            sxz00 = sxz00 + sxz * sxz;
            syz00 = syz00 + syz * syz;

            % Sample stress tensor correlation function
            tstress = tstress + 1;
            if mod(tstress, ITSTRESS0) == 0
                % New t=0
                tstress0 = tstress0 + 1;
                ttel = mod(tstress0 - 1, T0Max) + 1;
                tt0(ttel) = tstress;
                sxy0(ttel) = sxy;
                sxz0(ttel) = sxz;
                syz0(ttel) = syz;
            end
            for t = 1:min(tstress0, T0Max)
                dt = tstress - tt0(t) + 1;
                if dt < TMAx
                    nstress(dt) = nstress(dt) + 1;
                    sxyt(dt) = sxyt(dt) + sxy * sxy0(t);
                    sxzt(dt) = sxzt(dt) + sxz * sxz0(t);
                    syzt(dt) = syzt(dt) + syz * syz0(t);
                end
            end
        end

        if mod(NTVACF, Is) == 0
            tvacf = tvacf + 1;
            % Sample velocity auto-correlation function and mean square displacement

            if mod(IT0, tvacf) == 0
                % New t=0
                t0 = t0 + 1;
                ttel = mod(t0 - 1, T0Max) + 1;
                ttv0(ttel) = tvacf;
                for i = 1:NPART
                    rx0(i, ttel) = X(i);
                    ry0(i, ttel) = Y(i);
                    rz0(i, ttel) = Z(i);
                    vxt0(i, ttel) = VX(i);
                    vyt0(i, ttel) = VY(i);
                    vzt0(i, ttel) = VZ(i);
                end
            end
            for t = 1:min(t0, T0Max)
                dt = tvacf - ttv0(t) + 1;
                if dt < TMAx && dt * dtime <= TDIFMAX
                    display("ASIGNACION DE NVAC")
                    nvacf(dt) = nvacf(dt) + 1
                    r2asum = 0;
                    for i = 1:NPART
                        vacf(dt) = vacf(dt) + VX(i) * vxt0(i, t) + VY(i) * vyt0(i, t) + VZ(i) * vzt0(i, t);
                        r2a = (X(i) - rx0(i, t))^2 + (Y(i) - ry0(i, t))^2 + (Z(i) - rz0(i, t))^2;
                        r2t(dt) = r2t(dt) + r2a;
                        r2asum = r2asum + r2a;
                    end
                    r2asum = r2asum / NPART;
                    % Print mean square displacement to file for t=1,10,100,etc
                    nbl = 1;
                    iout = 49;
                    while nbl <= dt
                        iout = iout + 1;
                        if nbl - dt + 1 == 0
                            %fprintf(iout, '%f %f\n', (dt - 1) * dtime, r2asum);
                        end
                        nbl = 10 * nbl;
                    end
                end
            end
        end
    elseif Switch == 0
        % Initialize
        % Radial distribution function:
        ngr = 0;
        nhgr = NHIsmax;
        delg = HBOX / nhgr;
        delgi = 1 / delg;
        g = zeros(1, nhgr);
        
        % Stress tensor
        tstress0 = 0;
        dstresstime = NSAMP * Delt * IGR;
        sxy00 = 0;
        sxz00 = 0;
        syz00 = 0;
        
        % Velocity auto-correlation function and mean square displacement
        t0 = 0;
        tvacf = 0;
        dtime = NSAMP * Delt * NTVACF;
        r2t = zeros(1, TMAx);
        nvacf = zeros(1, TMAx);
        vacf = zeros(1, TMAx);
        nstress = zeros(1, TMAx);
        sxyt = zeros(1, TMAx);
        sxzt = zeros(1, TMAx);
        syzt = zeros(1, TMAx);

    elseif Switch == 2
        % Writes results to disk
        tempav = tempav / ngr;
        fac = 4 * PI ./ (3 * NPART * delg^3);
        %fprintf(69, '%d %f %f %f\n', NPART, tempav, rho, tempav);
        for i = 1:nhgr
            r = i * delg;
            rm = r - 0.5 * delg;
            r2 = r * r;
            g(i) = g(i) * fac / (r2 * r2);
            %fprintf(69, '%f %f\n', rm, g(i));
        end
        
        fac = 3 * BOX^3 / (2 * NSAMP * Delt * IGR);
        for t = 1:TMAx
            sxyt(t) = sxyt(t) * fac / (rho * temp * nstress(t));
            sxzt(t) = sxzt(t) * fac / (rho * temp * nstress(t));
            syzt(t) = syzt(t) * fac / (rho * temp * nstress(t));
            %fprintf(67, '%f %f %f %f\n', (t - 1) * dstresstime, sxyt(t), sxzt(t), syzt(t));
        end
        
        % Velocity auto-correlation function and mean square displacement
        fac = 1 / (NPART * dtime);
        for t = 1:TMAx
            vacf(t) = vacf(t) * fac / nvacf(t);
            r2t(t) = r2t(t) / nvacf(t);
            %display(68, '%f %f %f\n', (t - 1) * dtime, vacf(t), r2t(t));
        end
        r2tRet = r2t;
    end
end
