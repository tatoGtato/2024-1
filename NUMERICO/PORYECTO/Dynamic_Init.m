clear all 
clc 

molecularDynamics(T, dt, n, d)

function molecularDynamics(T, dt, n, d)

   [x,xm] = dynamic_init(n, temp, dt);
   t = 0;

   while (t < T)
       [f, en] = force(n, x, rc);
       integrate(f, en, n, x, xm);
       t = t + dt;
   end 
end 

function [x,xm] = dynamic_init(n, temp, dt)
    sumv = 0;
    sumv2 = 0;
    x = zeros(n);
    v = zeros(n);
    xm = zeros(n);
    
    for i = 1 : n
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

function [x, y] = lattice_pos(i, n)
    % Convert the linear index to 2D coordinates
    y = ceil(i / n);       % Row index
    x = mod(i-1, n) + 1;   % Column index
end

function [f, en] = force(n, x, rc)
    en = 0;
    f = zeros(n);

    for i = 1 : n-1
        for j = i+1 : n
            xr = x(i) - x(j);
            xr = xr-box*nint (xr/box);
            r2 = xr^2;

            if (r2 < rc*2)
                r2i = 1/r2;
                r6i = r2i^3;
                ff = 48 * r2i * r6i * (r6i - 0.5);
                f(i) = f(i) + ff*xr;
                f(j) = f(j) - ff*xr;
                en = en + 4 * r6i * (r6i -1) - 4((1/rc^12) - (1/rc^6))
            end
        end
    end
end

function [etot, temp] = integrate(f, en, n, x, xm)
    sumv = 0;
    sumv2 = 0;

    for i = 1 :n 
        xx = 2*x(i) - xm(i) + delt^2 * f(i);
        vi = (xx-xm(i))/(2*delt);
        sumv = sumv + vi;
        sumv2 =sumv2 + vi^2;
        xm(i) = x(i);
        x(i) = xx;
    end 
    temp = sumv2/(3*n);
    etot = (en + 0.5*sumv2)/n;
end 




