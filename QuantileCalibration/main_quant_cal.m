%% Define the interpolating part here
lower_middle_quantile = .2;
upper_middle_quantile = .8;

%% Define some probabilities
p_lower = .01:.01:lower_middle_quantile;
p_middle = lower_middle_quantile:.01:upper_middle_quantile;
p_upper = upper_middle_quantile:.01:.98;
p = [p_lower, p_middle, p_upper];

%% Get the quantiles with those probabilities + add noise
% rng(2372110);

q_lower  = -log(1-p_lower);
% q_lower  = norminv(p_lower);
q_lower  = q_lower + randn(1,length(q_lower)) * 0.02;

q_middle = -log(1-p_middle);
% q_middle = norminv(p_middle);
q_middle = q_lower(end) - q_middle(1) + 0.1 + q_middle + randn(1,length(q_middle)) * 0.03;

q_upper = -log(1-p_upper);
% q_upper = norminv(p_upper);
q_upper = q_middle(end) - q_upper(1) + 0.1 +q_upper + randn(1,length(q_upper)) * 0.05;
q = [q_lower, q_middle, q_upper];

%% Perform quantile calibration
quantile_calibration(p, q, .15, .85, true);
