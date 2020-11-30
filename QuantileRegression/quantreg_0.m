function [p, fval, tau] = quantreg_0(x, y, tau, noPrint, constraints)
% Fit best quantile regression model: q_y = alpha * x + x ^ beta * z(tau),
% where
% p(1) = alpha
% p(2) = beta
% p(2 + i) = z(tau(i))

%% Check input arguments
if nargin == 1 && length(x) == 4 && iscell(x)
    tau = x{3};
    y = x{2};
    x = x{1};
elseif nargin == 2
    if length(y) == 3
        constraints = y{3};
    else
        constraints = {};
    end
    tau = y{1};
    noPrint = y{2};
    y = x(:,2);
    x = x(:,1);
elseif nargin < 3
    error('Not enough input arguments.');
elseif nargin < 4
    noPrint = false;
    constraints = {};
elseif nargin < 5
    constraints = {};
end

if (any(tau<=0)) || (any(tau>=1))
    error('the percentile (tau) must be between 0 and 1.')
end

if size(x,1) ~= size(y,1)
    error('length of x and y must be the same.');
end

if numel(y) ~= size(y,1)
    error('y must be a column vector.')
end

if numel(x) ~= size(x,1)
    error('x must be a column vector.')
end

%% Get the Koenker rho function
flag_HT_median_constraints = false;
switch length(constraints)
    case 1
        if strcmp(constraints,'HT')
            optimising_rho = @(prm, x, y, tau, n_out, x_beta)opt_rho_HT(prm, x, y, tau, n_out, x_beta);
        end
    case 2
        if all(strcmp(constraints,'HT') + strcmp(constraints,'median'))
            optimising_rho = @(prm, x, y, tau, n_out, x_beta)opt_rho_HT_median(prm, x, y, tau, n_out, x_beta);
            flag_HT_median_constraints = true;
        end
end

%% Solve

tol = 1e-3;
phi = (1+sqrt(5))/2;
if ~flag_HT_median_constraints
    beta_hat = GoldenSectionSearch_BETA(@(b)evaluate_beta(b, x, y, tau, optimising_rho, 1), -1/phi, 1, noPrint, tol);
    % beta_hat = GridSearch_BETA(@(b)evaluate_beta(b, x, y, tau, optimising_rho, 1), -1, 1, noPrint, tol); toc();
    output = evaluate_beta(beta_hat, x, y, tau, optimising_rho, 2);
    
    alpha_hat = output{2};
    z_hat = output{3};
    
    p = [alpha_hat; beta_hat; z_hat];
    fval = output{1};
else
    beta_hat = GoldenSectionSearch_BETA(@(b)evaluate_beta_HT_median(b, x, y, tau, optimising_rho, 1), -1/phi, 1, noPrint, tol);
    % beta_hat = GridSearch_BETA(@(b)evaluate_beta(b, x, y, tau, optimising_rho, 1), -1, 1, noPrint, tol); toc();
    output = evaluate_beta_HT_median(beta_hat, x, y, tau, optimising_rho, 2);
    
    alpha_hat = output{2};
    c_hat = output{3};
    z_hat = output{4};
    
    p = [alpha_hat; c_hat; beta_hat; z_hat];
    fval = output{1};
end

end
