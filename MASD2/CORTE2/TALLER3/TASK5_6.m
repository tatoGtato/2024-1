% % Task 5
% 
% a = 0;  % Lower bound of the interval
% b = 1;   % Upper bound of the interval
% T = 1;     % Time
% N = 1e3;   % Number of samples
% 
% plot_histogram(a, b, T, N);
% 
% a = 0;  % Lower bound of the interval
% b = 1;   % Upper bound of the interval
% T = 1;     % Time
% N = 1e3;   % Number of samples
% 
% probability = compute_probability(a, b, T, N);
% disp(['Probability of finding a point inside the interval (a, b) at T = 1: ', num2str(probability)]);
% 
% function plot_histogram(a, b, T, N)
%     % Function to plot a histogram of random points
%     points = random_walk(T, N);
% 
%     % Plot histogram
%     histogram(points, 'Normalization', 'probability');
%     hold on;
% 
%     % Plot vertical lines indicating interval (a, b)
%     line([a, a], ylim, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
%     line([b, b], ylim, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
% 
%     xlabel('Position');
%     ylabel('Probability');
%     title('Histogram of Random Points');
%     legend('Histogram', 'Interval (a, b)');
%     hold off;
% end
% 
% function probability = compute_probability(a, b, T, N)
%     % Function to compute the probability of finding a point inside the interval (a, b) at T = 1
%     points = random_walk(T, N);
%     count = sum(points > a & points < b);  % Count points inside (a, b)
%     probability = count / N;
% end
% 
% function points = random_walk(T, N)
%     % Function to perform a random walk
%     x = rand(1, N) - 0.5;  % Initial positions
%     for t = 1:100*T  % Discretizing time
%         x = x + sqrt(2*T/100) * randn(1, N);  % Random walk
%     end
%     points = x;
% end


%%
Task 6
function plot_histogram(T, N)
    % Function to perform Monte Carlo simulation and plot histogram of final positions
    num_steps = 100 * T;
    final_positions = zeros(1, N);

    for i = 1:N
        % Perform random walk for each particle
        x = rand - 0.5;  % Initial position
        for t = 1:num_steps
            x = x + sqrt(2*T/100) * randn;  % Random walk
        end
        final_positions(i) = x;  % Store final position
    end

    % Plot histogram of final positions
    histogram(final_positions, 'Normalization', 'pdf');
    xlabel('Position');
    ylabel('Probability Density');
    title(['Final Positions (T = ', num2str(T), ')']);

    % Plot theoretical Gaussian distribution
    hold on;
    x_values = linspace(min(final_positions), max(final_positions), 1000);
    pdf_values = normpdf(x_values, mean(final_positions), std(final_positions));
    plot(x_values, pdf_values, 'r', 'LineWidth', 2);
    legend('Monte Carlo Simulation', 'Gaussian Distribution');
    hold off;
end

T = 1;     % Time
N = 1e5;   % Number of particles

plot_histogram(T, N);



