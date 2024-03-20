function [Pnew,i] = newtonRaphson(f, fp, p0, TOL, No)
    i = 1;
    p = p0;
    while i <= No
        Pnew = p - (f(p)/fp(p));
        if (abs(Pnew - p) < TOL)
            return
        end
        i = i+1;
        p = Pnew;
    end
end
