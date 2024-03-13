%%
%Aproximada

Uo = zeros(1,20);


K = 1/2;
dx = pi/20;
dt = pi/20;
muuuu = (K*dt)/(dx^2);

%Construimos Uo
i = 1;
for x = 0.0: dx :pi
    if (x < pi/2)
        Uo(i) = x;

    elseif (x == pi/2)
        Uo(i) = pi/2;

    else
        Uo(i) = pi-x;
    end
    i = i +1;
end


%Construimos la aproximacion
for j = 1: 1 : 21
    U = Uo.*(1 - 4*(muuuu)*sin((K*dt)/2)^2);
    Uo = U;
end


%%
%EXACTA
bk = @(k) (4*(-1)^((k+1)/2))/pi*k^2;

Dx = pi/20;

Dt = pi/20;

U = @(x,t, k) bk(k)*sin(k*x)*exp(-k^2*t)