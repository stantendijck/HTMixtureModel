%% QuantileCalibration

% Get z and tau by doing some calculations
thr = 0.7;
tau = [0.05:0.05:0.7,0.71:0.01:0.95]';

org_tau_high = tau(tau > thr);
mrg_tau_high = (org_tau_high - min(org_tau_high)) / (1-min(org_tau_high));

eps = 1e-6;
mrg_tau_high = mrg_tau_high + eps;

mu = 5; sigma = 1; xi = 0.1;

z_high = mu + sigma * ((1-mrg_tau_high).^(-xi)-1)/xi;
z_high = z_high + randn(size(z_high))*0.03;

for i = 2:length(z_high)
    if z_high(i) < z_high(i-1)
        z_high(i) = z_high(i-1);
    end
end

org_tau_low = tau(tau <= thr);
z_low = norminv(org_tau_low)*3 + 3.3;
z_low = z_low + randn(size(z_low)) * 0.06;

for i = 2:length(z_low)
    if z_low(i) < z_low(i-1)
        z_low(i) = z_low(i-1);
    end
end

%% Perform the quantile calibration
z = [z_low; z_high];

output_tau = (0.05:0.001:0.95)';

n = 100;
xi_estimates = zeros(n,1);
thr = linspace(0.6,0.9,n);

tic();
for ithr = 1:n
    [~,opt_p] = quantile_calibration(z, tau, thr(ithr), output_tau, ithr == 1);
    xi_estimates(ithr) = opt_p(3);
end
toc();

figure(4); clf;
plot(thr,xi_estimates)






