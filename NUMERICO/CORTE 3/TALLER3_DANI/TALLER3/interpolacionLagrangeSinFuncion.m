%recordar y1 = linspace(-5,5,7) rango -5 a 5 7 puntos incluyendo -5 y 5
%para evaluar
%para la interpolacion inversa igualar a 0 utilizar metodo anterior
function c = interpolacionLagrangeSinFuncion(puntosX, puntosY)
    syms x
    polinomioLagrange = 0;
    n = size(puntosX, 2);
    
    for i = 1:n
        numerador = 1;
        denominador = 1;
        
        for j = 1:n  
            if i ~= j
                numerador = numerador * (x - puntosX(j));
                denominador = denominador * (puntosX(i) - puntosX(j));
            end
            if i == j
                yAMultiplicar = puntosY(j);
            end
        end
        
        polinomioLagrange = polinomioLagrange + yAMultiplicar*(numerador / denominador);
    end
    
    c = simplify(polinomioLagrange, 'Steps', 50);
    
    
end