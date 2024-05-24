function coeficientes = ajustePolinomial(X, Y, grado)
    A = zeros(length(X), grado + 1);
    for i = 1:grado + 1
        A(:, i) = X.^(grado + 1 - i);
    end

    % Resolver el sistema de ecuaciones normales para obtener los coeficientes
    coeficientes = (A' * A) \ (A' * Y');
end
