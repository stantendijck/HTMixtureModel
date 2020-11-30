function NLOGL = HT_mixture_likelihood(prm, x, y, rsdLike, logx)
% Calculate the complete Heffernan-Tawn likelihood
if nargin < 4
	logx = log(x);
end

%% Defaults
alpha = prm(1,:);
beta = prm(2,:);
mu = prm(3,:);
sig = prm(4,:);
xbeta = exp(logx .* beta);

noc = length(alpha);

if any(sig < 0) || any(abs(alpha) > 1) || any(beta > 1)
	NLOGL = inf;
	return
end

%% KeefConstraints
kcond = ones(2,1);
if 1
    kcond(1) = HeffernanTawn.KeefConstraints(x, y, [alpha;beta], 1, xbeta);
	kcond(2) = HeffernanTawn.KeefConstraints(x, y, [alpha;beta], 0, xbeta);
end

if ~all(kcond)
	NLOGL = Inf;
	return
end

%% Get the allocation probabilities
allocprob = zeros(size(rsdLike));
for i = 1:noc
	allocprob(:,i) = rsdLike(:,i) ./ sum(rsdLike,2);
end

%% Calculate the likelihood for all x separately
L = sum(allocprob ./ (sqrt(2*pi) * sig .* xbeta) .* exp(-(y - alpha .* x - mu .* xbeta).^2 ./ (2 * (sig .* xbeta).^2)),2);

%% Combine the likelihoods
NLOGL = - sum(log(L));


end