a = 0;

b = 2;

i = 0;

fa = f2(a);

tol = 10^(-2);

maxIteraciones = 10;


while i <= maxIteraciones
    p = (a+b)/2;
    fp = f2(p);
    if (fp == 0 || (b-a)/2 < tol)
            p
            break
    end 
    i = i+1;
    
    if (fa*fp > 0)
        a = p;
        fa = fp;
    else 
        b = p;
    end 
end



function f1 = f1(x)
    f1 = x^3 + 4*x^2 - 10;
end

function f2 = f2(x)
    f2 = exp(x) - 4 + x;
end
