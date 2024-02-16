x0 = 2;
x1 = 2.75;
x2 = 4;

puntos  = [2 2.75 4]
fPuntos  = [0 -3 1]
iteracionesMax = length(puntos)


for i = 1:iteracionesMax
    L = 1 

    for j= 1:iteracionesMax
        if j ~= i 
            L = L * (x_eval - x(j)) / (x(i) - x(j))
        end
    end

    y_eval = y_eval + y(i) * L 
    
end 
