function [p,i] = biseccion(f, a, b, tol, No)
    i = 0;
    fa = f(a);
    
    while i <= No
        p = (a+b)/2;
        fp = f(p);
        if (fp == 0 || (b-a)/2 < tol)
                return
                
        end 
        i = i+1;
        
        if (fa*fp > 0)
            a = p;
            fa = fp;
        else 
            b = p;
        end 
    end
end



