clc

% Ejemplo de uso
x = 1:10;
f = log(x);
x_estimado = 1.3;


lagrangeInterpolation(x,f,x_estimado)


% function valorIntermedio = lagrangeInterpolation(x, y, x_eval)
%   % Convertir a valores de coma flotante
%   x = double(x);
%   y = double(y);
%   x_eval = double(x_eval);
% 
%   n = length(x); % Número de puntos
% 
%   % Detectar si x_eval está fuera del rango
%   if x_eval < x(1) || x_eval > x(n)
%     disp('Error: El valor a estimar está fuera del rango de datos.');
%     return;
%   end
% 
%   % Interpolación cuadrática para el primer intervalo
%   if x_eval <= x(2)
%     i = 1;
%     L = ones(1, n);
%     for j = 1:n
%       if j ~= i
%         L(j) = L(j) * (x_eval - x(j)) / (x(i) - x(j));
%       end
%     end
%     valorIntermedio = y(i) * L(i);
%     return;
%   end
% 
%   % Interpolación cúbica para los intervalos intermedios
%   if x_eval > x(2) && x_eval < x(n-1)
%     for i = 2:n-1
%       if x_eval >= x(i) && x_eval < x(i+1)
%         break;
%       end
%     end
%     L = ones(1, n);
%     for j = 1:n
%       if j ~= i && j ~= i+1
%         L(j) = L(j) * (x_eval - x(j)) / (x(i) - x(j));
%       end
%     end
%     valorIntermedio = y(i) * L(i) + y(i+1) * L(i+1) * (x_eval - x(i)) / (x(i+1) - x(i));
%     return;
%   end
% 
%   % Interpolación cuadrática para el último intervalo
%   if x_eval >= x(n-1)
%     i = n;
%     L = ones(1, n);
%     for j = 1:n
%       if j ~= i
%         L(j) = L(j) * (x_eval - x(j)) / (x(i) - x(j));
%       end
%     end
%     valorIntermedio = y(i) * L(i);
%   end
% end

function valorIntermedio = lagrangeInterpolation(x, y, x_eval)
  % Convertir a valores de coma flotante
  x = double(x);
  y = double(y);
  x_eval = double(x_eval);

  n = length(x); % Número de puntos

  % Detectar si x_eval está fuera del rango
  if x_eval < x(1) || x_eval > x(n)
    disp('Error: El valor a estimar está fuera del rango de datos.');
    return;
  end

  % Interpolación cúbica para todos los puntos
  L = ones(1, n); % Vector de polinomios base de Lagrange
  for j = 1:n
        for k = 1:n
          if k ~= j
            L(j) = L(j) * (x_eval - x(k)) / (x(j) - x(k));
          end
        end
  end

  valorIntermedio = 0;
  for i = 1:n
    valorIntermedio = valorIntermedio + y(i) * L(i);
  end
end

