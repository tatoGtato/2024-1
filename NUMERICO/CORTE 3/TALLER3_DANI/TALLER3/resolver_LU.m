function x = resolver_LU(L, U, b)
    % Inicializar el vector solución y con ceros del mismo tamaño que b
    y = zeros(size(b));
    
    % Obtener el número de elementos en b (debe ser igual a las filas de L y U)
    n = length(b);
    
    % Primera fase: resolver Ly = b para y
    for i = 1:n
        % Calcular el i-ésimo elemento de y utilizando la fila i de L
        % y los elementos previamente calculados de y
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end
    
    % Inicializar el vector solución x con ceros del mismo tamaño que b
    x = zeros(size(b));
    
    % Segunda fase: resolver Ux = y para x
    for i = n:-1:1
        % Calcular el i-ésimo elemento de x utilizando la fila i de U
        % y los elementos previamente calculados de x
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end