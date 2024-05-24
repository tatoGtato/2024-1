function c = coeficienteCorrelacion(puntosX, puntosY, puntosYCalculados)
% Definición de la función coeficienteCorrelacion
% Calcula el coeficiente de correlación entre los valores observados y los valores calculados.

    sr = 0; % Inicializa la suma de los cuadrados de las diferencias entre los valores observados y calculados
    st = 0; % Inicializa la suma de los cuadrados de las diferencias entre los valores observados y la media de los valores observados

    mediaYi = mean(puntosY); % Calcula la media de los valores observados

    % Bucle para calcular las sumas de los cuadrados de las diferencias
    for i=1:length(puntosX)
        sr = sr + (puntosYCalculados(i) - puntosY(i))^2; % Suma de los cuadrados de las diferencias entre los valores observados y calculados
        st = st + (puntosY(i) - mediaYi)^2; % Suma de los cuadrados de las diferencias entre los valores observados y la media de los valores observados
    end

    % Calcula el coeficiente de correlación utilizando la fórmula estándar
    c = sqrt((st - sr) / st); 

end
