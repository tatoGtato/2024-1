
p0 = 0.2;
p1 = 0.9;
tol = 10^(-5);
iteracionesMax = 20;

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
    p0 = p1;
    q0 = q1;
    p1 = p;
    q1 = f(p);
end

function f = f(x)
    % f = x^5-2*x^3 - log(x);
    f = x^3 + 2*x^2 + 10*x - 20;
end
