%x = @(z) (1+z^(-1)) / (z^-1 + 4*z^-2 + 0.5*z^-3);

num = [0 0 1 1];
den = [1.5 4 1 0];

[b,a,n,m] = eqtflength(num,den)

[z,p,k] = tf2zpk(b,a)

zplane(b,a)
text(real(z)-0.1,imag(z)-0.1,"Zeros")
text(real(p)-0.1,imag(p)-0.1,"Poles")


h = impz(b,a,10)