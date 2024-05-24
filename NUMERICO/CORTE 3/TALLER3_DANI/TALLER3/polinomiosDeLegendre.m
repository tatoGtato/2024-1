function c = polinomiosDeLegendre(f, polinomio, numeroPolinomioLegendre)
    %RECORDAR QUE ESTO FUNCIONA EN EL INTERVALO -1 1
    syms x 

    if polinomio == "legendre"
        % Inicializa un arreglo de celdas para almacenar los polinomios de Legendre
        polinomiosLegendre = {};
        polinomiosLegendre{1} = 0; % P0(x) = 0 (aunque típicamente P0(x) = 1)
        polinomiosLegendre{2} = 1; % P1(x) = x

        % Genera los polinomios de Legendre hasta el número especificado
        for i = 0:numeroPolinomioLegendre-1
            % Usa la recurrencia de los polinomios de Legendre para generar el siguiente polinomio
            polinomiosLegendre{end+1} = (1/(i+1))*((2*i+1)*x*polinomiosLegendre{end} - i * polinomiosLegendre{end-1});
        end

        % Elimina el primer elemento (P0(x) = 0) para que el primer polinomio sea P1(x) = x
        polinomiosLegendre = polinomiosLegendre(2:end);

        % Inicializa un arreglo de celdas para almacenar los coeficientes de los polinomios
        aj = {};
        tamano = size(polinomiosLegendre);
        iteracion = tamano(2);

        % Ciclo para calcular los coeficientes de los polinomios de Legendre
        for i = 0:iteracion-1
            % Calcula la norma del polinomio de Legendre
            norma = 2/(2*i+1);

            % Obtiene el polinomio de Legendre actual
            g = polinomiosLegendre{i+1};

            % Define la función a integrar como el producto de 'f' y el polinomio de Legendre
            funcionAIntegrar = f*g;

            % Convierte la función simbólica a una función anónima de MATLAB para integración numérica
            ht = matlabFunction(funcionAIntegrar);

            % Parámetros de integración
            a = -1;
            b = 1;
            n = 10;
            result = 0;
            paso = (b-a)/n;

            % Aplica la regla de Simpson para la integración numérica
            for j = a:paso:b
                x0 = j;
                x1 = j + paso/2;
                x2 = j + paso;

                % Ajusta el último punto si es necesario para no exceder el intervalo
                if x2 > b
                    x2 = b;
                end

                % Acumula el resultado de la integración
                result = result + integracionReglaSimpson(x0, x1, x2, ht);
            end

            % Calcula el coeficiente 'aj' y lo almacena en el arreglo
            aj{end+1} = (1/norma)*result;
        end

        % Construye el polinomio de aproximación como combinación lineal de los polinomios de Legendre
        polinomio = 0;
        for i = 1:iteracion
            polinomio = polinomio + polinomiosLegendre{i} * aj{i};
        end

        % Convierte el resultado a una aproximación numérica con 'vpa' (Variable Precision Arithmetic)
        c = vpa(polinomio);
    end
end