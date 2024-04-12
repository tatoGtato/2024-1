%%
%a
clc
clear all

y = @(x) tan(x^3);
x = 3;
h = 0.1;

"Punto medio 3"
derivadaPuntoMedioTresPuntos(y, x, h)

"Punto estremo 3"
derivadaExtremoTresPuntos(y, x, h)

"Punto medio 5"
derivadaCincoPuntosPuntoMedio(y, x, h)

"Punto estremo 5"
derivadaCincoPuntosExtremo(y, x, h)

%%
%b
clc
clear all

y = @(x) exp(x) + x;
x = 2;
h = 0.2;

"Punto medio 3"
derivadaPuntoMedioTresPuntos(y, x, h)

"Punto estremo 3"
derivadaExtremoTresPuntos(y, x, h)

"Punto medio 5"
derivadaCincoPuntosPuntoMedio(y, x, h)

"Punto estremo 5"
derivadaCincoPuntosExtremo(y, x, h)

%%
%3 puntos medio

function derivada = derivadaPuntoMedioTresPuntos(func, x, h)
    % func: la función de la cual se quiere calcular la derivada
    % x: el punto en el cual se quiere calcular la derivada
    % h: el tamaño de paso utilizado para la aproximación
    
    % Evaluar la función en los dos puntos necesarios
    % f_x_h = func(x + h);
    % f_x_mh =func(x - h);


    % Evaluar la función en los dos puntos necesarios con redondeo a 5
    % digitos ejemplo 3
    f_x_h = round(func(x + h), 5);
    f_x_mh = round(func(x - h), 5);

    % Aplicar la fórmula del punto medio de tres puntos
    % derivada = (f_x_h - f_x_mh) / (2*h);

    % Aplicar la fórmula del punto medio de tres puntos con redondeo
    % ejemplo 3
    derivada = round((f_x_h - f_x_mh) / (2*h), 5) ;
end


%%
%3 puntos extremo
function derivada = derivadaExtremoTresPuntos(func, x, h)
    % func: la función de la cual se quiere calcular la derivada
    % x: el punto en el cual se quiere calcular la derivada
    % h: el tamaño de paso utilizado para la aproximación

    % Evaluar la función en los tres puntos necesarios
    f_x = func(x);
    f_x_h = func(x + h);
    f_x_2h = func(x + 2*h);

    % Aplicar la fórmula del punto extremo de tres puntos
    derivada = (-3*f_x + 4*f_x_h - f_x_2h) / (2*h);
end

%%
%5 puntos medio
function derivada = derivadaCincoPuntosPuntoMedio(func, x, h)
    % func: la función de la cual se quiere calcular la derivada
    % x: el punto en el cual se quiere calcular la derivada
    % h: el tamaño de paso utilizado para la aproximación

    % Evaluar la función en los cinco puntos necesarios
    f_xm2h = func(x - 2*h); % f(x - 2h)
    f_xmh = func(x - h);   % f(x - h)
    f_xph = func(x + h);   % f(x + h)
    f_xp2h = func(x + 2*h); % f(x + 2h)

    % Aplicar la fórmula de cinco puntos en el punto medio
    derivada = (f_xm2h - 8*f_xmh + 8*f_xph - f_xp2h) / (12*h);
end

%%
%5 puntos extremo
function derivada = derivadaCincoPuntosExtremo(func, x, h)
    % func: la función de la cual se quiere calcular la derivada
    % x: el punto en el cual se quiere calcular la derivada
    % h: el tamaño de paso utilizado para la aproximación

    % Evaluar la función en los cinco puntos necesarios
    f_x = func(x);        % f(x)
    f_xh = func(x + h);   % f(x + h)
    f_x2h = func(x + 2*h); % f(x + 2h)
    f_x3h = func(x + 3*h); % f(x + 3h)
    f_x4h = func(x + 4*h); % f(x + 4h)

    % Aplicar la fórmula de cinco puntos en el extremo
    derivada = (-25*f_x + 48*f_xh - 36*f_x2h + 16*f_x3h - 3*f_x4h) / (12*h);
end

