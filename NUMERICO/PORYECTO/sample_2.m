function r2t = sample_2(Switch, Is, En, Vir, Enk, Delt)
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

        if mod(Is, IGR) == 0
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

        if mod(Is, NTVACF) == 0
            tvacf = tvacf + 1;
            % Sample velocity auto-correlation function and mean square displacement
            if mod(tvacf, IT0) == 0
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
                    nvacf(dt) = nvacf(dt) + 1;
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
                            fprintf(iout, '%f %f\n', (dt - 1) * dtime, r2asum);
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
            fprintf(69, '%f %f\n', rm, g(i));
        end
        
        fac = 3 * BOX^3 / (2 * NSAMP * Delt * IGR);
        for t = 1:TMAx
            sxyt(t) = sxyt(t) * fac / (rho * temp * nstress(t));
            sxzt(t) = sxzt(t) * fac / (rho * temp * nstress(t));
            syzt(t) = syzt(t) * fac / (rho * temp * nstress(t));
            fprintf(67, '%f %f %f %f\n', (t - 1) * dstresstime, sxyt(t), sxzt(t), syzt(t));
        end
        
        % Velocity auto-correlation function and mean square displacement
        fac = 1 / (NPART * dtime);
        for t = 1:TMAx
            vacf(t) = vacf(t) * fac / nvacf(t);
            r2t(t) = r2t(t) / nvacf(t);
            display(68, '%f %f %f\n', (t - 1) * dtime, vacf(t), r2t(t));
        end
    end
end