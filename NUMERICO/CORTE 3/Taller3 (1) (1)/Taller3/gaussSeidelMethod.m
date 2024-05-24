function [x, success] = gaussSeidelMethod(n, A, b, xo, tol, N)
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
            % Inicializar las sumatorias
            sumatoria1 = 0;
            sumatoria2 = 0;
            
            % Calcular la primera sumatoria (para j < i)
            for j = 1:i-1
                sumatoria1 = sumatoria1 + A(i,j) * x(j);
            end
            
            % Calcular la segunda sumatoria (para j > i)
            for j = i+1:n
                sumatoria2 = sumatoria2 + A(i,j) * xo(j);
            end
            
            % Calcular el nuevo valor de x(i)
            x(i) = (b(i) - sumatoria1 - sumatoria2) / A(i,i);
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
