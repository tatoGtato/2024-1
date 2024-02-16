clc 

n = 200;
t = 1;
c = 0.5;


f = @(x) 0.4*exp(-300*(x-0.5)^2) + 0.1*exp(-300*(x-0.65)^2);
U0_X = zeros(n,1);
shiftMatrix = zeros(n,n);
row = 1;
I = eye(n);


%MATRIZ DE U CON T = 0 
for i = 0.0 :0.005: t
   if i ~= 1
    U0_X(row, 1) = f(i);
    row = row + 1; 
   end
end


%MATRIZ DE FORWARD
for i = 1.0: 1.0 :n
   if i ~= 200 
    shiftMatrix(i, i+1) = 1;
   end
end

%U = -c*(shiftMatrix*U0_X) + (c+1)*U0_X

for i = 0.0: 0.005 : t
    U = -c*(shiftMatrix*U0_X) + (c+1)*U0_X;
    U0_X = U;
end

U

% figure
% X = linspace(0,1,200);
% Y = 0.4*exp(-300*(X-0.5).^2) + 0.1*exp(-300*(X-0.65).^2);
% stem(X,Y)
% 
% 
% exp(1)