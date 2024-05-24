function [error_estimado, coef_correlacion] = calcularErrorYCorrelacion(X, Y, Y_aprox)
    % Calcular el error estimado (RMSE)
    error_estimado = sqrt(mean((Y - Y_aprox).^2));

    % Calcular el coeficiente de correlación (R)
    R = corrcoef(Y, Y_aprox);
    coef_correlacion = R(1, 2);
end
