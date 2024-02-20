%%PLOT
x = linspace(-10,10);
plot(x, exp(x).*cos(x) - x.^2 + 3*x);
ylim([-5 5])
axes0
