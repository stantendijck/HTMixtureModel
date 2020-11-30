function [lower_quantiles, upper_quantiles] = get_upper_and_lower_quantile_functions(QR, j_noc)

%% Better estimation procedures
rho = @(z, t)(z .* (t - (z < 0)));

upper_quantiles = cell(j_noc, 1);
lower_quantiles = cell(j_noc, 1);

cps = [0;QR.q{1,j_noc};1];
for ic = 1:j_noc
    curr_alpha = QR.p{1,j_noc}{ic}(1);
    curr_c = QR.p{1,j_noc}{ic}(2);
    curr_beta = QR.p{1,j_noc}{ic}(3);
    curr_x_beta = QR.x .^ curr_beta;
    
    tau_lower = cps(ic);
    tau_upper = cps(ic+1);
    
    fun = @(z, t)(sum(rho(QR.y - curr_alpha * QR.x - curr_c - z * curr_x_beta, t)));
    z_lower = fminsearch(@(z)(fun(z, tau_lower)), 0);
    z_upper = fminsearch(@(z)(fun(z, tau_upper)), 0);
    
    lower_quantiles{ic} = @(x)(curr_c + curr_alpha * x + z_lower * x .^ curr_beta);
    upper_quantiles{ic} = @(x)(curr_c + curr_alpha * x + z_upper * x .^ curr_beta);
end

end
