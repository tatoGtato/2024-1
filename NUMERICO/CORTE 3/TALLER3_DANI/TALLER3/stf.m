% La función stf (Short-Time Fourier Transform) calcula la serie de Fourier de una
% función f en un intervalo de tiempo T, considerando N términos de la serie y 
% límites inferior (lower_lim) y superior (upper_lim) para la integración. 
% La serie de Fourier es una representación de una función como una suma infinita 
% de funciones seno y coseno, y es útil para analizar señales periódicas.
function fourier_series = stf(f, T, N, lower_lim, upper_lim)
    L = T / 2; % Semiperíodo de la función
    omega = 2 * pi / T; % Frecuencia angular de la función

    % Ajusta los límites de integración si son -1, lo que indica los límites por defecto
    if lower_lim == -1
        lower_lim = -L;
    end
    if upper_lim == -1
        upper_lim = L;
    end
    
    % Calcula el coeficiente a0 de la serie de Fourier utilizando la regla del trapecio compuesta
    a0 = (1 / L) * ReglaTrapecioCompuesta(f, lower_lim, upper_lim, 1000);

    % Inicializa los vectores para los coeficientes an y bn
    an = zeros(1, N);
    bn = zeros(1, N);
    
    % Calcula los coeficientes an y bn utilizando la regla del trapecio compuesta
    for n = 1:N
        an(n) = (1 / L) * ReglaTrapecioCompuesta(@(x) f(x) .* cos(n * omega * x), lower_lim, upper_lim, 1000);
        bn(n) = (1 / L) * ReglaTrapecioCompuesta(@(x) f(x) .* sin(n * omega * x), lower_lim, upper_lim, 1000);
    end
    
    % Inicializa la serie de Fourier con la mitad del coeficiente a0
    fourier_series = @(x) a0 / 2;

    % Construye la serie de Fourier sumando los términos de coseno y seno ponderados por an y bn
    for n = 1:N
        fourier_series = @(x) fourier_series(x) + an(n) * cos(n * omega * x) + bn(n) * sin(n * omega * x);
    end
end