x1=1e-4:1e-4:10;
x2=3:1e-4:10;

for m=0:2
    Jm=besselj(m,x1);
    figure(1)
    plot(x1,Jm);
    hold on;
    Km=besselk(m,x2);
    figure(2)
    plot(x2,Km);
    hold on;
end
figure(1);
xlabel('x');
ylabel('J_m(x)');
legend('m=0','m=1','m=2','m=3','m=4','m=5');
hold off;
figure(2);
xlabel('x');
ylabel('K_m(x)');
legend('m=0','m=1','m=2','k=3','m=4','m=5');
hold off;