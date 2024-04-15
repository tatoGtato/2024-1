
function dydt = ODE3a(t, y, k)
    dydt = -k * sqrt(y);
end