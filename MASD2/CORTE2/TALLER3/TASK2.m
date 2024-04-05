%%
%Aproximada
clc
clear all 

K = 1;
dx = pi/20;
dt = (pi^2/800);
muuuu = (K*dt)/(dx^2);

T = (3*pi^2/80);
n = 20+1;



Uo = zeros(1,n);
shiftMatrixF = zeros(n,n);
shiftMatrixB = zeros(n,n);


%MATRIZ DE FORWARD
for i = 1.0: 1.0 :n
   if i ~= n 
    shiftMatrixF(i, i+1) = 1;
   end
end


%MATRIZ DE BACKWARD
for i = 1.0: 1.0 :n
   if i ~= n 
    shiftMatrixB(i+1, i) = 1;
   end
end


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
    i = i + 1;
end

%Construimos la aproximacion

for j = 0: dt : T
    U = (muuuu) * (Uo * shiftMatrixF) + (1 - 2*muuuu)*Uo + muuuu * muuuu * (Uo * shiftMatrixB);
    Uo = U;
end

figure
X = linspace(0, pi, n);
Y = Uo;

size(X)
size(Y)

stem(X,Y)

%% 
clc
clear all 

K = 1;
dx = pi/20;
dt = (pi^2/800);
mu = (K*dt)/(dx^2);

T = (3*pi^2/80);
n = 20+1;

Uo = zeros(1,n);

% Construct initial condition Uo
for i = 1:n
    x = (i-1) * dx;
    if (x < pi/2)
        Uo(i) = x;
    elseif (x == pi/2)
        Uo(i) = pi/2;
    else
        Uo(i) = pi-x;
    end
end

% Construct the approximation using central difference method
A = zeros(n,n);
for i = 2:n-1
    A(i,i-1) = mu;
    A(i,i) = 1 - 2*mu;
    A(i,i+1) = mu;
end
A(1,1) = 1;
A(n,n) = 1;

% Construct the next state U using implicit method
for j = 0:dt:0.2
    U = A * Uo';
    Uo = U';
end

% Plotting
figure
X = linspace(0, pi, n);
Y = Uo;
stem(X,Y)



%%
%EXACTA
% Constants
L = pi;
T = (3 * pi^2) / 80;
dx = pi / 20;
dt = pi^2 / 1000;

% Number of points
Nx = round(L / dx);
Nt = round(T / dt);

% Initialize solution matrix
u = zeros(Nx, Nt);

% Function to compute bk
compute_bk = @(k) 4 * (-1)^((k+1)/2) / (pi * k^2);

% Compute the solution
x = linspace(0, L, Nx);
t = linspace(0, T, Nt);

for i = 1:Nt
    for j = 1:Nx
        for k = 1:100  % summing up to a large enough value
            u(j, i) = u(j, i) + compute_bk(2*k - 1) * sin((2*k - 1) * x(j)) * exp(-((2*k - 1)^2) * t(i));
        end
    end
end
u = -u; 
% Plotting
figure;
for i = 1:5:Nt  % Plot at different times
    plot(x, u(:, i), 'DisplayName', ['t = ', num2str(t(i))]);
    hold on;
end

title('Numerical Solution of the PDE at Different Times');
xlabel('x');
ylabel('u(x, t)');
legend('Location', 'best');
grid on;
hold off;
