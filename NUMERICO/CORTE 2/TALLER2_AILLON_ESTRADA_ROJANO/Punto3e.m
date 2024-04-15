%%
%A
clc
clear all

f = @(x,y) 4*exp(0.8*x)-0.5*y;
a = 0;
b = 2;
h = 0.2;
yo = 2;
e = 10^(-5);

hmin = 0.1 * h;
hmax = 0.4;

[T, W, H] = rungeKuttaFehlberg(f, a, b, yo, e, hmax, hmin)

options = odeset('RelTol', e);
[tOde, wOde] = ode45(f, [a b], yo, options);

tiledlayout(2,1)
% Top plot
nexttile
plot(T, W);
title('Aproximacion')
xlim([0 2])

% Bottom plot
nexttile
plot(tOde,wOde)
title('Exacta')
xlim([0 2])

%%
%B
clc
clear all

f = @(x,y) 10*exp( -((x-2)^2) / (2*(0.075)^2) );
a = 0;
b = 4;
h = 0.25;
yo = 0.5;
e = 10^(-5);

hmin = 0.1 * h;
hmax = 0.75;

[T, W, H] = rungeKuttaFehlberg(f, a, b, yo, e, hmax, hmin)

options = odeset('RelTol', e);
[tOde, wOde] = ode45(f, [a b], yo, options);

tiledlayout(2,1)
% Top plot
nexttile
plot(T, W);
title('Aproximacion')
xlim([0 4])

% Bottom plot
nexttile
plot(tOde,wOde)
title('Exacta')
xlim([0 4])


%%
function [T, W, H] = rungeKuttaFehlberg(f, a, b, alpha, TOL, hmax, hmin)
    % Estimación inicial del número máximo de pasos
    maxSteps = ceil((b - a) / hmin) + 1;
    
    t = a;
    w = alpha;
    h = hmax;
    FLAG = 1;
    stepCount = 1;

    % Preasignación con la estimación del tamaño máximo
    T = zeros(1, maxSteps);
    W = zeros(1, maxSteps);
    H = zeros(1, maxSteps);
    
    T(stepCount) = t;
    W(stepCount) = w;
    H(stepCount) = h;

    while FLAG == 1
        K1 = h * f(t, w);
        K2 = h * f(t + 1/4*h, w + 1/4*K1);
        K3 = h * f(t + 3/8*h, w + 3/32*K1 + 9/32*K2);
        K4 = h * f(t + 12/13*h, w + 1932/2197*K1 - 7200/2197*K2 + 7296/2197*K3);
        K5 = h * f(t + h, w + 439/216*K1 - 8*K2 + 3680/513*K3 - 845/4104*K4);
        K6 = h * f(t + 1/2*h, w - 8/27*K1 + 2*K2 - 3544/2565*K3 + 1859/4104*K4 - 11/40*K5);

        R = abs(1/h * (1/360*K1 - 128/4275*K3 - 2197/75240*K4 + 1/50*K5 + 2/55*K6));

        if R <= TOL
            t = t + h;
            w = w + 25/216*K1 + 1408/2565*K3 + 2197/4104*K4 - 1/5*K5;
            
            stepCount = stepCount + 1;
            T(stepCount) = t;
            W(stepCount) = w;
            H(stepCount) = h;
        end

        delta = 0.84 * (TOL/R)^(1/4);
        if delta <= 0.1
            h = 0.1 * h;
        elseif delta >= 4
            h = 4 * h;
        else
            h = delta * h;
        end

        if h > hmax
            h = hmax;
        end

        if t >= b
            FLAG = 0;
        elseif t + h > b
            h = b - t;
        elseif h < hmin
            FLAG = 0;
            fprintf('Minimum h exceeded\n');
        end
    end

    % Recorta los vectores T, W, H al tamaño real utilizado
    T = T(1:stepCount);
    W = W(1:stepCount);
    H = H(1:stepCount);
end