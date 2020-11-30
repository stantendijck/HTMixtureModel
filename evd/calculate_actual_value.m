function actual_value = calculate_actual_value(model, method, xl, xu, yl, yu)

switch method
    case 'simulation'
        [large_x, large_y] = DataSml(model.name, model.large_n, model.tau_cp, false);
        actual_value = sum(xl < large_x & large_x < xu & yl < large_y & large_y < yu)/model.large_n * 1/2 * exp(-min(large_x));
    case 'math'
        switch model.name
            case 'log'
                cdf_alog_gumbel_margins = @(x, alpha)( ...
                    exp(-(exp(-x(1)/alpha) + exp(-x(2)/alpha)).^alpha));
                alog_prm.alpha = 0.5;
                
                cdf = @(x)(cdf_alog_gumbel_margins(Gumbel_iCDF(Laplace_CDF(x)), alog_prm.alpha));
                
                actual_value = cdf([xu,yu]) - cdf([xl,yu]) - cdf([xu,yl]) + cdf([xl,yl]);
            case 'alog'
                cdf_alog_gumbel_margins = ...
                    @(x, theta, alpha)( ...
                exp(-theta(1) * exp(-x(1)) - theta(2) * exp(-x(2)) - ...
                    ((1-theta(1))^(1/alpha) * exp(-x(1)/alpha) + ...
                     (1-theta(2))^(1/alpha) * exp(-x(2)/alpha)).^alpha));
                % gumbel_icdf = @(u)(-log(-log(u)));
                % laplace_cdf = @(x)(1/2 + 1/2 * sign(x) .* (1 - exp(-abs(x))));
                alog_prm.theta = [model.tau_cp, 0.5];
                alog_prm.alpha = 0.5;
                
                cdf = @(x)(cdf_alog_gumbel_margins(Gumbel_iCDF(Laplace_CDF(x)), alog_prm.theta, alog_prm.alpha));
                
                actual_value = cdf([xu,yu]) - cdf([xl,yu]) - cdf([xu,yl]) + cdf([xl,yl]);
            case 'normal'
                rho = 0.5;
                
                Mu = zeros(2, 1);
                Sig = eye(2);
                Sig(1, 2) = rho;
                Sig(2, 1) = rho;
                
                xl_n = norminv(Laplace_CDF(xl));
                xu_n = norminv(Laplace_CDF(xu));
                yl_n = norminv(Laplace_CDF(yl));
                yu_n = norminv(Laplace_CDF(yu));                
               
                actual_value = mvncdf([xl_n; yl_n], [xu_n; yu_n], Mu, Sig);
            case 'normal_log_mixture'
                rho = 0.5;
                alpha = 0.5;
				prob = model.tau_cp;
               
                actual_value = MVNLOG_CDF([xu;yu], prob, rho, alpha) - MVNLOG_CDF([xu;yl], prob, rho, alpha) - MVNLOG_CDF([xl;yu], prob, rho, alpha) + MVNLOG_CDF([xl;yl], prob, rho, alpha);
            otherwise
                error("calculate_actual_value(model, method) only supported for model.name = 'alog' if method == 'math'");
        end
    otherwise
        error("calculate_actual_value(model, method) only supported for method = 'math' or 'simulation'");
end
