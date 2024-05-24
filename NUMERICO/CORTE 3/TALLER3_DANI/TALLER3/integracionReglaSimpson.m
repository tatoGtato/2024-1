function c = integracionReglaSimpson(x0, x1, x2, f)

    h = (x2-x0)/2;
    
    c = (h/3)*(f(x0)+4*f(x1)+f(x2));
end