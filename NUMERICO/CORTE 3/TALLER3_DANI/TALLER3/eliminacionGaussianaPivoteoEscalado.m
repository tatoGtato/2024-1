%está estructurado como una función que toma una matriz matriz y un vector z, y devuelve el vector solución x.
function x = eliminacionGaussianaPivoteoEscalado(matriz, z)
    % Obtiene el tamaño de la matriz de coeficientes y el vector de resultados
    tamano = size(matriz);
    tamanoZ = size(z);
    
    % Calcula los factores de escalamiento como el máximo valor absoluto por fila
    scalingFactors = max(abs(matriz), [], 2); 
    
    % Proceso de eliminación gaussiana con pivoteo escalado
    for i = 1:tamano(2)
        % Encuentra el valor máximo escalado en la columna actual y su fila pivote
        [maxValue, pivotRow] = max(abs(matriz(i:end, i)) ./ scalingFactors(i:end));
        pivotRow = pivotRow + i - 1;  % Ajusta el índice de la fila pivote
        
        % Intercambia la fila actual con la fila pivote si es necesario
        if pivotRow ~= i
            matriz([i, pivotRow], :) = matriz([pivotRow, i], :);
            z([i, pivotRow]) = z([pivotRow, i]);
        end
        
        % Realiza la eliminación para hacer ceros debajo de la diagonal principal
        for j = (i + 1):tamano(1)
            factor = matriz(j, i) / matriz(i, i);  % Calcula el factor de eliminación
            matriz(j, :) = matriz(j, :) - factor * matriz(i, :);  % Actualiza la fila j
            z(j) = z(j) - factor * z(i);  % Actualiza el vector de resultados
        end
    end
    
    % Inicializa el vector de solución x con ceros
    x = zeros(1, tamanoZ(1));
    
    % Proceso de sustitución hacia atrás para obtener la solución del sistema
    for i = tamanoZ(1):-1:1
        x(i) = z(i);  % Inicializa x(i) con el valor correspondiente en z
        for j = (i + 1):tamanoZ(1)
            x(i) = x(i) - matriz(i, j) * x(j);  % Resta los términos conocidos
        end
        x(i) = x(i) / matriz(i, i);  % Divide por el coeficiente diagonal para obtener x(i)
    end
end