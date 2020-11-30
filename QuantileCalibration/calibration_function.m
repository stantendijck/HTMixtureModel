function y = calibration_function(Prm, input_p, input_q, uniform_p)

if Prm(1) < 0
    y = 1e8;
    return;
end

opt_mu = min(input_q);

gpd_cdf = @(x, sigma, xi)((1 - (1 + xi * (x - opt_mu)/sigma) .^ (-1 / xi)) .* (x < opt_mu - sigma/xi | xi > 0) + 1e9 * (xi < 0 & x >= opt_mu - sigma/xi));
% inverse_gpd_cdf = @(u, sigma, xi)(opt_mu + sigma/xi * ((1-u) .^(-xi) - 1));

extended_input_p = [input_p(1) - (input_p(2) - input_p(1)),input_p,input_p(end) + (input_p(end) - input_p(end-1))];

eps = 1e-8;
y = 0;
for i = 1:length(input_p)
    d_input_p = (extended_input_p(i+2) - extended_input_p(i))/2;
    
    p_model = gpd_cdf(input_q(i) + eps, Prm(1), Prm(2));
    y = y + d_input_p * (p_model - uniform_p(i))^2 / (input_p(i) * (1 - input_p(i)));
end

%% Horizontal difference
% y = sum((inverse_gpd_cdf(uniform_p, Prm(1), Prm(2)) - input_q).^2);

%% Vertical difference
% eps = 1e-8;
% y = sum((gpd_cdf(input_q + eps, Prm(1), Prm(2)) - uniform_p).^2);

end