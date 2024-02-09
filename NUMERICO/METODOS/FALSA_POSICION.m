clc

p0 = 1;
p1 = 2;
tol = 10^(-5);
iteracionesMax = 500;

i = 2;
q0 = f(p0);
q1 = f(p1);

while i <= iteracionesMax
    p = p1 - ( (q1*(p1-p0))/(q1-q0) );
    if abs(p - p1) < tol 
        p
        break
    end
    i = i+1;
    q = f(p);
    if (q * q1 < 0)
        p0 = p1;
        q0 = q1;
    end
    p1 = p;
    q1 = q;
    
end




function f = f(x)
    %f = x^3 + 4*x^2 - 10;
    f = x^10 - 10;
end
