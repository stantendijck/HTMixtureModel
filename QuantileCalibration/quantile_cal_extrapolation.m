function output_q = quantile_cal_extrapolation(input_q, input_p, output_p, tp)

if isempty(input_p)
    output_q = [];
    return;
end

if nargin < 4
    if abs(input_p(1)) > abs(1 - input_p(end))
        tp = 'upper';
    else
        tp = 'lower';
    end
end

eps = 1e-4;
switch tp
    case 'upper'
        range_p = 1 - input_p(1);
        uniform_p = (input_p - input_p(1))/range_p;
        
        output_range_p = 1 - output_p(1);
        output_uniform_p = (output_p - output_p(1))/output_range_p;
    case 'lower'
        input_p = 1 - input_p(end:-1:1);
        range_p = 1 - input_p(1);
        uniform_p = (input_p - input_p(1))/range_p;
        uniform_p = uniform_p + eps;
        
        input_q = -input_q(end:-1:1);
        
        output_p = 1 - output_p(end:-1:1);
        output_range_p = 1 - output_p(1);
        output_uniform_p = (output_p - output_p(1))/output_range_p;
end

output_uniform_p = min(1-1e-5,max(1e-5, output_uniform_p));
            

opt_mu = min(input_q);

inverse_gpd_cdf = @(u, sigma, xi)(opt_mu + sigma/xi * ((1-u) .^(-xi) - 1));

calibration_function2 = @(Prm)(calibration_function(Prm, input_p, input_q, uniform_p));
[opt_prm, ~, exitflag, output] = fminsearch(calibration_function2, [1, 1]); %,optimset('FunValCheck','on','Display','iter-detailed'));
if exitflag == 0
    debug = true;
end

gpd_cdf2 = @(x, sigma, xi)((1 - (1 + xi * (x - opt_mu)/sigma) .^ (-1 / xi)) .* (x < opt_mu - sigma/xi | xi > 0) + 1 * (xi < 0 & x >= opt_mu - sigma/xi));

output_q = inverse_gpd_cdf(output_uniform_p, opt_prm(1), opt_prm(2));

% figure; plot(input_p, input_q);
% hold on;
% plot(output_p, output_q);

if strcmp(tp,'lower')
    output_q = -output_q(end:-1:1);
end

end