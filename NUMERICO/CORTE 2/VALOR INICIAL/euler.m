b = 0;
a = 2;
alpha = 0.5;
N = 10;
yp = @(y,t) y - t^2 + 1;

[t,w] = eulers(a,b,alpha,N,yp)

function [t,w] = eulers(a,b,alpha,N,f)
    h = (b-a)/N;
    t = a;
    w = alpha;
    for i = 1.0: 1 :N
        w = w + h*f(t,w);
        t = a + h;
        disp("EN " + i)
        disp("t " + t)
        disp("w " + w)
    end
    return 
end