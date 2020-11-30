function v = MVNLOG_iCDF(p, prob, rho, alpha)
% Calculate the inverse cdf of a mixture of a bivariate normal and logistic on Laplace margins such that
% P(X > v & Y > v) = p
% or 
% P(X > v(1) & Y < v(2)) = p(1); P(X > v(1)) = p(2);

% Multivariate Gaussian parameters
Mu = zeros(2, 1);
Sig = eye(2);
Sig(1, 2) = rho;
Sig(2, 1) = rho;

if length(p) == 1
	cdf = @(v)(MVNLOG_CDF([v; v], prob, rho, alpha));
	fun = @(v)(1 - 2*Laplace_CDF(v) + cdf(v));
	
	v = fminsearch(@(v)(abs(log(fun(v)) - log(p))),10);
else
	% Calculate r
	r = Laplace_iCDF(1-p(1)/p(2));
	
	% Calculate v0
	cdf = @(v)(MVNLOG_CDF([r; v], prob, rho, alpha));
	fun = @(v)(log(Laplace_CDF(v) - cdf(v)));
	v0 = fminsearch(@(w)(abs(fun(w) - log(p(1)))),5);
	
	v = [r, v0];
end






end