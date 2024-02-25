%% PUNTO 3
clc 

n = 200;
t = 1;
c = 0.5;

%PUNTO 3
f = @(x) 0.4*exp(-300*(x-0.5).^2) + 0.1*exp(-300*(x-0.65).^2);

U0_X = zeros(n,1);
shiftMatrix = zeros(n,n);
row = 1;
I = eye(n);


%MATRIZ DE U CON T = 0 
for i = 0.0 :0.005: t
   if i ~= t
    U0_X(row, 1) = f(i);
    row = row + 1; 
   end
end


%MATRIZ DE FORWARD
for i = 1.0: 1.0 :n
   if i ~= n 
    shiftMatrix(i, i+1) = 1;
   end
end

%U = -c*(shiftMatrix*U0_X) + (c+1)*U0_X
%ONE SIDED METHOD
for i = 0.0: 0.005 : t
    U = -c*(shiftMatrix*U0_X) + (c+1)*U0_X;
    U0_X = U;
end

figure
X = linspace(0,1,200);
if (t == 0)
    Y = f(X);
else 
    Y = abs(U);
end
stem(X,Y)


%% PUNTO 5

n = 200;
t = 0.3;
c = 0.5;

%PUNTO 3
f = @(x) 0.4*exp(-300*(x-0.5).^2) + 0.1*exp(-300*(x-0.65).^2);

%PUNTO 6    

% f = @(x) 1 - (x-0.7)/(0.1) %|x-0.7| <= 0.1
% 

U0_X = zeros(n,1);
U0_Xm1 = zeros(n,1);
shiftMatrix = zeros(n,n);
row = 1;
I = eye(n);


%MATRIZ DE U CON T = 0 
for i = 0.0 :0.005: t
   if i ~= t
    U0_X(row, 1) = f(i);
    row = row + 1; 
   end
end


%MATRIZ DE FORWARD
for i = 1.0: 1.0 :n
   if i ~= n 
    shiftMatrix(i, i+1) = 1;
   end
end

%U = -c*(shiftMatrix*U0_X) + (c+1)*U0_X
% LAX-WENDROFF METHOD
for i = 0.0: 0.005 : t
    U = 1/2*c*(c-1)*(shiftMatrix*U0_X) - (c^2 - 1)*U0_X + 1/2 * c*(c+1) * U0_Xm1; 
    U0_X = U;
    U0_Xm1 = U0_X;
end


% figure
% X = linspace(0,1,200);
% Y = f(X);
% stem(X,Y)

figure
X = linspace(0,1,200);
if (t == 0)
    Y = f(X);
else 
    Y = abs(U);
end
stem(X,Y)

%% PUNTO 6

n = 200;
t = 0.3;
c = 0.5;

 f = @(x) 1 - (x-0.7)/(0.1) %|x-0.7| <= 0.1


U0_X = zeros(n,1);
U0_Xm1 = zeros(n,1);
shiftMatrix = zeros(n,n);
row = 1;
I = eye(n);


%MATRIZ DE U CON T = 0 
for i = 0.0 :0.005: t
   if i ~= t
    U0_X(row, 1) = f(i);
    row = row + 1; 
   end
end


%MATRIZ DE FORWARD
for i = 1.0: 1.0 :n
   if i ~= n 
    shiftMatrix(i, i+1) = 1;
   end
end

%U = -c*(shiftMatrix*U0_X) + (c+1)*U0_X
% LAX-WENDROFF METHOD
for i = 0.0: 0.005 : t
    U = 1/2*c*(c-1)*(shiftMatrix*U0_X) - (c^2 - 1)*U0_X + 1/2 * c*(c+1) * U0_Xm1; 
    U0_X = U;
    U0_Xm1 = U0_X;
end


% figure
% X = linspace(0,1,200);
% Y = f(X);
% stem(X,Y)

figure
X = linspace(0,1,200);
if (t == 0)
    Y = f(X);
else 
    Y = abs(U);
end
stem(X,Y)