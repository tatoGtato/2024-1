function coeficientes = ajusteLegendre(X, Y, n)
    syms x
    l = lagrangeCoefficients(X,Y);
    Xn= (x+15/9).*(9/2); % para x in [-1,1]
    f = simplify(subs(l,x,Xn));

    x0 = -1; xf = 1;

    phi{1} = 0;
    phi{2} = 1;
    for i = 0:n-1
        phi{end+1} = (1/(i+1))*((2*i+1)*x*phi{end} - i * phi{end-1});
    end
    phi = phi(2:end);

    for j = 0:n
        aux= int(f.* phi(j+1));
        in = vpa(subs(aux,xf)-subs(aux,x0),5);
        aj(j+1) = vpa((2*j + 1)* in/2);
    end

    coeficientes = aj;
end


