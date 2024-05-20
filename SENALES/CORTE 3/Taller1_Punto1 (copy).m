clear all
clc

f1=100;
f2=200;
fs= 1000;
t=0:1/fs: 0.049;

%SEÑAL
x = 2*cos(2*pi*f1*t) + 1.5*sin(2*pi*f2*t);

N = size(x,2);
k = 0:1/(N-1):1;

%ARMONICOS
Ck = getCk(N, x);

%PASA BAJA
pb = Pb(0.4);
%PASA ALTA
pa = Pa(0.1);
%BANDA ANCHA
ba = Ba(0.1,0.1);
%RESONADOR
Fr = fr(0.49,0.49);

%SEÑALES FILTRADAS
xFiltradaBaja = Ck.*pb;
xFiltradaAlta = Ck.*pa;

nexttile
stem(t,x);
title("Señal")

nexttile
stem(k,abs(Ck));
title("Armonicos")

nexttile
plot((0:100)./100.*0.5,abs(pb))
title("Pasa baja")

nexttile
plot((0:100)./100.*0.5,abs(pa))
title("Pasa alta")

nexttile
plot((0:100)./100.*0.5,abs(pa))
title("Banda Ancha")

nexttile
plot((0:100)./100.*0.5,abs(Fr))
title("Resonador")

nexttile
plot(k,xFiltradaBaja);
title("Señal filtrada PB")

nexttile
plot(k,xFiltradaAlta);
title("Señal filtrada PA")

%%
%Fucniones Usadas

function Hpb = Pb(fc)
        n=-25:25;
        h = 2*fc.*sinc(2*fc.*n);
        m=(0:100)./100;
        Hpb =freqz(h,1,pi.*m);
    return
end

function Hpa = Pa(fc)
        n=-25:25;
        h = (-1).^n*2*fc.*sinc(2*fc.*n);
        m=(0:100)./100;
        Hpa =freqz(h,1,pi.*m);
    return
end

function Hba = Ba(fc1, fc2)
        n=-25:25;
        h = (-1).^n*2*fc2.*sinc(2*fc2.*n) - 2*fc1.*sinc(2*fc1.*n);
        m=(0:100)./100;
        Hba =freqz(h,1,pi.*m);
    return
end

function Fr = fr(r,w)
        h = r*exp(-1i*w);
        m=(0:100)./100;
        Fr =freqz(h,1,pi.*m);
    return
end

function Fr = fr(r,w)
        h = r*exp(-1i*w);
        m=(0:100)./100;
        Fr =freqz(h,1,pi.*m);
    return
end

function Ck = getCk(N, x)
    Ck = zeros(N,1);
    for k = 1 : N
        Ck(k) = (1/N)*sum(x.*exp(-1i*2*pi*(k-1).*(0:N-1)./N));
    end 
    return 
end 

function Xn = getXn(N, Ck)
    Xn = zeros(N,1); 
    for k = 1 : N
        if (k < ceil(N/2)) 
            k
            Xn(k) = sum(Ck'.*exp(-1i*2*pi*(k-1).*(0:N-1)./N ));
        else
            Xn(k) = sum(Ck'.*exp(-1i*2*pi*(k-1).*(0:N-1)./N ));
        end
    end 
    Xn(1)
    return 
end 