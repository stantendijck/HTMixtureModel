function Y = opt_rho_HT(prm, x, y, tau, n_out, x_beta)
% prm = [alpha,beta], we penalise abs(alpha) < 1; beta <= 1 is implied by
% implementation method

tol = 1e-3;

alpha = prm(1);
beta = prm(2);

if nargin < 5
    x_beta = x.^beta;
end

factor = 1 + 100*(~all(HeffernanTawn.KeefConstraints(x,y,[alpha;beta],[0,1],x_beta)));

opt_z = zeros(length(tau),1);
fval = zeros(length(tau),1);

z_start = -10; z_end = 10;
tau_start = tau(1); tau_end = tau(end);

fun = @(z)(sum(abs((y - (alpha*x + z*x_beta)) .* (tau_start - (y - (alpha*x + z*x_beta)<0)))));
opt_z(1) = GoldenSectionSearch_Z(fun, z_start, z_end, tol);

if length(tau) > 1
    fun = @(z)(sum(abs((y - (alpha*x + z*x_beta)) .* (tau_end - (y - (alpha*x + z*x_beta)<0)))));
    opt_z(end) = GoldenSectionSearch_Z(fun, opt_z(1), opt_z(1) + z_end - z_start, tol);
end

if length(tau) > 2
    for i = 2:length(tau)-1
        fun = @(z)(sum(abs((y - (alpha*x + z*x_beta)) .* (tau(i) - (y - (alpha*x + z*x_beta)<0)))));
        opt_z(i) = GoldenSectionSearch_Z(fun, opt_z(i-1), opt_z(end), tol);
    end
end

for i = 1:length(tau)
    fun = @(z)(sum(abs((y - (alpha*x + z*x_beta)) .* (tau(i) - (y - (alpha*x + z*x_beta)<0)))));
    fval(i) = fun(opt_z(i));
end

if length(tau) <= 2 || n_out == 2
    w = ones(1,length(tau));
else
    begin_tau = max(0,3/2*tau(1) - tau(2)/2);
    end_tau = min(1,3/2*tau(end) - tau(end-1)/2);
    w = diff([begin_tau,(tau(2:end) + tau(1:end-1))/2,end_tau])/min(diff(tau));
end

if n_out == 1
    Y = factor * (w * fval); % for optimising it
else
    Y = {factor * (w * fval), opt_z};   % for getting the actual value.
end

end