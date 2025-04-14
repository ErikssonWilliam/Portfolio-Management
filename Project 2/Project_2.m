%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = readmatrix('portfolioManagerV4Project 2 2025.xls', 'Sheet', 'assetStat', 'Range', 'B3:B12');
Cov = readmatrix('portfolioManagerV4Project 2 2025.xls', 'Sheet', 'assetStat', 'Range', 'S4:AB13');
r = readmatrix('portfolioManagerV4Project 2 2025.xls', 'Sheet', 'Refinitiv', 'Range', 'K5:K5');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_mean_rate = mu - r;
C_inverse_diff_Mean_rate = Cov\diff_mean_rate; %Take inverse and mulitply

one_vector = ones(length(mu), 1);
% Tangency portfolio
denominator_M = one_vector' * C_inverse_diff_Mean_rate;
xM = C_inverse_diff_Mean_rate / denominator_M;

%Min-var portfolio (x_R)
C_inverse_one = Cov \ one_vector;
denominator_R = one_vector' * C_inverse_one;
xR = C_inverse_one / denominator_R;

%x_T
Q_matrix = inv(Cov) - (C_inverse_one * C_inverse_one') / denominator_R;
xT = xR + Q_matrix * mu;

% Tangency portfolio metrics
VaRxM = xM' * Cov * xM; % Variance (scalar)
SigmaxM = sqrt(VaRxM);
MuxM = mu' * xM; % Expected return (scalar)

% Beta calculation 
beta = Cov * xM / SigmaxM; % (10x1)

dec = 4;

CheckResults(dec, Cov, mu, r, xR, xT, xM, MuxM, SigmaxM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Efficient Frontier, Analytic approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_values = [linspace(-20,-0.5,100), linspace(0.5,20,100)];
num_gamma = length(gamma_values);

x_gamma = zeros(length(mu), num_gamma);
mu_bar = zeros(1, num_gamma);
sigma_squared = zeros(1, num_gamma);


for i = 1:num_gamma
    if gamma == 0
    continue; % skip or handle specially
    end
    gamma = gamma_values(i);
    % Compute x_gamma
    x_gamma(:, i) = (1/gamma) * xT + (1 - 1/gamma) * xR;
    
    % Expected return
    mu_bar(i) = mu' * x_gamma(:, i);
    
    % Variance (more efficient calculation)
    sigma_squared(i) = x_gamma(:, i)' * Cov * x_gamma(:, i);
end

sigma = sqrt(sigma_squared);

denominator_sigma_alt = sqrt(diff_mean_rate' * (Cov \ diff_mean_rate));
sigma_alt = abs(mu_bar - r) / denominator_sigma_alt;

figure;
hold on;
plot(sigma, mu_bar, 'b-', 'LineWidth', 1.5); 
plot(sigma_alt, mu_bar, 'g-', 'LineWidth', 1.5);
scatter(SigmaxM, MuxM, 100, 'r', 'filled'); % Tangency portfolio point
xlabel('Standard Deviation');
ylabel('Expected Return ');
title('Efficient Frontier: Expected return vs. Standard Deviation');
grid on;
legend('Efficient Frontier','Efficient Frontier (Alt)', 'Location', 'best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Efficient Frontier, Numeric approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First optimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimization settings
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');

% Equality constraint (sum of weights = 1)
Aeq = one_vector';  % Budget constraint
beq = [1];

% No inequality constraints
A = [];
b = [];

% Initial guess (equal-weighted portfolio)
x0 = one_vector/length(mu);

% Define gamma values (risk aversion parameters)
gamma_values = linspace(0.5, 20, 100);
num_gamma = length(gamma_values);

% Preallocate results
x_num = zeros(length(mu), num_gamma);  % Optimal weights
mu_num = zeros(1, num_gamma);          % Expected returns
sigma_num = zeros(1, num_gamma);       % Standard deviations
min_var_portfolio = zeros(1, num_gamma);       

for i = 1:num_gamma
    
    % Solve optimization problem
    [x_opt, obj] = fmincon(@(x) portfolio_objective_function(x, mu, gamma_values(i), Cov), x0, A, b, Aeq, beq, [], [], [], options);
    
    % Store results
    x_num(:, i) = x_opt;
    mu_num(i) = mu'*x_opt;
    min_var_portfolio(i) = portfolio_minvar(x_opt,Cov);
    sigma_num(i) = sqrt(min_var_portfolio(i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second optimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Optimization settings
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');

% Equality constraint (sum of weights = 1)
Aeq = [one_vector',1];  % Budget constraint, z = 1'*x + y
beq = [1];

% No inequality constraints
A = [];
b = [];

% Initial guess (equal-weighted portfolio)
x0 = [one_vector/length(mu); 0];

% Preallocate results
x_num2 = zeros(length(mu), num_gamma);  % Optimal weights
mu_num2 = zeros(1, num_gamma);          % Expected returns
sigma_num2 = zeros(1, num_gamma);       % Standard deviations
min_var_portfolio2 = zeros(1, num_gamma);       

for i = 1:num_gamma
    
    % Solve optimization problem
    [z_opt, obj] = fmincon(@(z) portfolio_objective_function2(r, z, mu, gamma_values(i), Cov), x0, A, b, Aeq, beq, [], [], [], options);
    x_opt = z_opt(1:end-1);
    % Store results
    x_num2(:, i) = x_opt;
    mu_num2(i) = mu'*x_opt;
    min_var_portfolio2(i) = portfolio_minvar(x_opt,Cov);
    sigma_num2(i) = sqrt(min_var_portfolio2(i));
end

figure;
hold on;

% First optimization frontier
plot(sigma_num, mu_num, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Numerical Frontier 1');

% Second optimization frontier
plot(sigma_num2, mu_num2, 'b-.', 'LineWidth', 1.5, 'DisplayName', 'Numerical Frontier 2');

% Mark tangency portfolio
scatter(SigmaxM, MuxM, 100, 'r', 'filled', 'DisplayName', 'Tangency Portfolio');

xlabel('Standard Deviation (\sigma)');
ylabel('Expected Return (\mu)');
title('Efficient Frontier Comparison: Optimization 1 vs 2');
grid on;


function f = portfolio_objective_function(x, mu, gamma, Cov)
    f = -mu(:)'*x(:) + (gamma/2)*(x(:)'*Cov*x(:));
end

function f = portfolio_minvar(x,C)
f = x(:)'*C*x(:);
end

function f = portfolio_objective_function2(r, z, mu, gamma, Cov)
    x = z(1:end-1);
    y = z(end);
    f = -mu(:)'*x(:) + (gamma/2)*(x(:)'*Cov*x(:)) - y*r;
end