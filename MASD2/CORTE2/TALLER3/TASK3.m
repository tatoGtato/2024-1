clear all
clc 

f = @(x) (1/sigma * sqrt(2*pi))* exp(-1/2 * ((x - miu)/sigma )^2)

a = -0.5; % Lower bound of the interval
b = 0.5;  % Upper bound of the interval
N = 10000; % Number of steps
dt = 0.01; % Time step

% Calculate the standard deviation of X
std_X = sqrt(N * dt);

% Calculate the probability using the CDF of the standard normal distribution
prob = normcdf(b / std_X) - normcdf(a / std_X);

disp(['Probability that the particle lies in the interval (', num2str(a), ', ', num2str(b), ') at T = 1: ', num2str(prob)]);