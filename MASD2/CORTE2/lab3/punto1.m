%%
% Aliasing
clc
clear all

n = 2^3;
f = zeros(n,n);

for i = 1:n
    f(1,i) = 1;
    f(i,1) = 1;
end 
f

for k =1:n-1
    for j = 1:n-1
        f(k,j) = exp(-2*pi*1i/n)^j;
    end 
end 

f

% Extract the real and imaginary parts of the fourth column
fourth_column_real = real(f(:, 4));
fourth_column_imag = imag(f(:, 4));

% Extract the real and imaginary parts of the (N-3)th column
N_minus_3_column_real = real(f(:, n-3));
N_minus_3_column_imag = imag(f(:, n-3));

% Plot the real parts
subplot(2, 1, 1);
plot(1:n, fourth_column_real, 'b', 'LineWidth', 2);
hold on;
plot(1:n, N_minus_3_column_real, 'r--', 'LineWidth', 2);
xlabel('Index');
ylabel('Real Part');
title('Comparison of Real Parts');
legend('Fourth Column', 'N-3 Column');

% Plot the imaginary parts
subplot(2, 1, 2);
plot(1:n, fourth_column_imag, 'b', 'LineWidth', 2);
hold on;
plot(1:n, N_minus_3_column_imag, 'r--', 'LineWidth', 2);
xlabel('Index');
ylabel('Imaginary Part');
title('Comparison of Imaginary Parts');
legend('Fourth Column', 'N-3 Column');

% Adjust subplot spacing
sgtitle('Comparison of Fourth Column with (N-3) Column');