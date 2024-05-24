% realiza una aproximación trigonométrica discreta a un conjunto de datos
% x_data: un vector que contiene los valores de la variable independiente (x) de los datos.
% 
% y_data: un vector que contiene los valores de la variable dependiente (y) de los datos.
% 
% m: el número de puntos de datos. Aunque no se utiliza directamente en el código, se asume que m es la longitud de x_data e y_data.
% 
% n: el número de términos trigonométricos a considerar en la aproximación.


% La función calcula los coeficientes de Fourier (ak y bk) utilizando las fórmulas
% de la serie de Fourier para datos discretos. Luego, construye la aproximación trigonométrica Sn
% como una suma de funciones seno y coseno ponderadas por estos coeficientes.
% 
% Después de calcular la aproximación, la función crea una figura y traza los datos originales 
% y la aproximación trigonométrica en el mismo gráfico para comparar visualmente la calidad de la aproximación. 
% Se añaden una leyenda, un título y etiquetas para los ejes, y se activa la cuadrícula para mejorar la legibilidad del gráfico.
% 
% Finalmente, la función devuelve la aproximación trigonométrica Sn como salida.

function Sn = trigonometric_approx(x_data, y_data, m, n)
    % ak y bk son vectores para almacenar los coeficientes de Fourier
    ak = zeros(1, n+1);
    bk = zeros(1, n);

    % Cálculo de los coeficientes ak para k desde 0 hasta n
    for k = 0:n
        ak(k+1) = (1/m) * sum(y_data .* cos(k*x_data));
    end
    
    % Cálculo de los coeficientes bk para k desde 1 hasta n-1
    for k = 1:n-1
        bk(k) = (1/m) * sum(y_data .* sin(k*x_data));
    end
    
    % Inicialización de la aproximación trigonométrica Sn con los primeros términos
    Sn = ak(1)/2 + ak(n+1)*cos(n*x_data);
    
    % Suma de los términos restantes a la aproximación Sn
    for k = 1:n-1
        Sn = Sn + ak(k+1)*cos(k*x_data) + bk(k)*sin(k*x_data);
    end
    
    % Creación de una figura y gráfico de los datos originales y la aproximación
    figure;
    plot(x_data, y_data, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Datos Originales');
    hold on; 
    plot(x_data, Sn, '--', 'LineWidth', 2, 'DisplayName', 'Aproximación Trigonométrica');
    legend('Location', 'best'); 
    title(['Aproximación Trigonométrica Discreta con ' num2str(n) ' términos']); 
    xlabel('x'); 
    ylabel('y'); 
    grid on; 
    hold off; 
end

% La función trigonometric_approx que has definido anteriormente calcula los 
% coeficientes de Fourier y construye una aproximación trigonométrica de los 
% datos de entrada utilizando una suma finita de funciones seno y coseno.
% 
% En el código que has proporcionado, estás llamando a la función trigonometric_approx con los siguientes argumentos:
% 
% x_data: un vector con 64 puntos de datos equiespaciados en el intervalo de [0, 6.3].
% 
% y_data: un vector con los valores correspondientes a x_data. Estos valores parecen seguir un patrón periódico.
% 
% m: la longitud de x_data dividida por 2, lo que sugiere que m es el número de puntos de datos en medio ciclo de la función periódica.
% 
% n: el número de términos trigonométricos a considerar en la aproximación, que has establecido en 25.
% 
% Al llamar a trigonometric_approx(x_data, y_data, m, n), estás pidiendo a la función que calcule una aproximación 
% trigonométrica de los datos utilizando 25 términos de la serie de Fourier. La función devolverá la aproximación Sn 
% y también generará un gráfico que muestra los datos originales y la aproximación trigonométrica para comparar visualmente la calidad de la aproximación.