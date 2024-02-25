clc

f=1/50;
N=200;
n=0:N-1;
truncs = 3
res = 2/truncs;
tiledlayout(4,1)
err = 0:N-1;

%LISTA DE MULTIPLOS
mults=1:truncs;
for i = 1.0:1:truncs+1 
    
end


x = sin(2*pi*f*n);
truncX = 0:N-1;
roundX = 0:N-1;

% %cuantificación de la señal
%%TRUNCAMIENTO
for i = 1.0:1:N
    truncatedPoint = mod(x(i),res);
    err(i) = truncatedPoint;
    truncX(i) = x(i) - truncatedPoint;
end

%%REDONDEO
for i = 1.0:1:N
    truncatedPoint = mod(x(i),res);
    if (truncatedPoint == 0)
        roundX(i) == x(i)
    else if(truncatedPoint > )
    end

end


nexttile
stem(n,x,'r');
grid on
title("Señal discreta x(n)")

nexttile
stem(n,truncX,'r');
title("Señal truncada")

nexttile
stem(n,roundX,'r');
title("Señal roudneada")

nexttile
stem(n,err,'r');
title("Error")

