function v = MVN_iCDF(p, rho)
% Calculate the inverse cdf of a bivariate normal such that
% P(X > v & Y > v) = p
% or 
% P(X > v(1) & Y < v(2)) = p(1); P(X > v(1)) = p(2);

% Multivariate Gaussian parameters
Mu = zeros(2, 1);
Sig = eye(2);
Sig(1, 2) = rho;
Sig(2, 1) = rho;

if length(p) == 1
	fun = @(w)(log(mvncdf([w; w], [inf; inf], Mu, Sig)));
	v_gaussian_scale = fminsearch(@(w)(abs(fun(w) - log(p))),0);
	
	% Check if error is small enough
	[p0, err] = mvncdf([v_gaussian_scale; v_gaussian_scale], [inf; inf], Mu, Sig);
else
	% Calculate r
	r = - norminv(p(1)/p(2));
	
	% Calculate v0
	fun = @(w)(log(mvncdf([r; -inf], [inf; w], Mu, Sig)));
	v0 = fminsearch(@(w)(abs(fun(w) - log(p(1)))),0);
	
	v_gaussian_scale = [r, v0];
	
	% Check if error is small enough
	[p0, err] = mvncdf([v_gaussian_scale(1); -inf], [inf; v_gaussian_scale(2)], Mu, Sig);
end

v = Laplace_iCDF(normcdf(v_gaussian_scale));




end