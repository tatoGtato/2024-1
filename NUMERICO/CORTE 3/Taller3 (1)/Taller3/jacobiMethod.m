function [x, success] = jacobiMethod(n, A, b, xo, tol, N)
    % Inicializar k
    k = 1;
    
    % Inicializar xo_ant para comparar con xo
    xo_ant = xo;
    
    % Inicializar la solución x
    x = xo;
    
    % Iterar hasta alcanzar el número máximo de iteraciones
    while k <= N
        % Iterar sobre todas las ecuaciones
        for i = 1:n
            % Inicializar la suma de los términos excepto el término i-ésimo
            sumatoria = 0;
            
            % Calcular la suma de los términos excepto el término i-ésimo
            for j = 1:n
                if j ~= i
                    sumatoria = sumatoria + A(i,j) * xo(j);
                end
            end
            
            % Calcular el nuevo valor de x(i)
            x(i) = (b(i) - sumatoria) / A(i,i);
        end
        
        % Comprobar si se alcanzó la tolerancia
        if norm(x - xo_ant) < tol
            success = true;
            return;
        end
        
        % Actualizar xo_ant para la siguiente iteración
        xo_ant = x;
        
        % Incrementar k
        k = k + 1;
    end
    
    % Si se excedió el número máximo de iteraciones
    success = false;
end
