function p = MVNLOG_CDF(x, prob, rho, alpha)
% Evaluate the CDF of a mixture of a Normal(rho) (p = prob) and Logistic(alpha) (p = 1 - prob) at [x(1), x(2)] on Laplace margins

%% Multivariate Gaussian parameters
Mu = zeros(2, 1);
Sig = eye(2);
Sig(1,2) = rho;
Sig(2,1) = rho;

%% Transform Laplace scale to Frechet scale
z = Frechet_iCDF(Laplace_CDF(x));

%% Transform Frechet scale to Gaussian scale
if prob ~= 0
	xu_normal = norminv(Frechet_CDF(z(1)/prob));
	yu_normal = norminv(Frechet_CDF(z(2)/prob));
else
	xu_normal = inf;
	yu_normal = inf;
end

%% Transform Frechet scale to Gumbel scale
if prob ~= 1
	xu_log = Gumbel_iCDF(Frechet_CDF(z(1)/(1-prob)));
	yu_log = Gumbel_iCDF(Frechet_CDF(z(2)/(1-prob)));
else
	xu_log = inf;
	yu_log = inf;
end

%% Calculate the Gaussian part
p1 = mvncdf([xu_normal;yu_normal],Mu,Sig);

%% Calculate the logistic part
logcdf = @(x, alp)(exp(-(exp(-x(1)/alp) + exp(-x(2)/alp)).^alp));
p2 = logcdf([xu_log;yu_log],alpha);

%% Combine the results
p = p1 * p2;


end