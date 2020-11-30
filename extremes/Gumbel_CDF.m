function u = Gumbel_CDF(x, mu, beta)
%Compute the cumulative distribution function of a Gumbel random variable
%with parameters (mu, beta) -> default = (mu = 0, beta = 1)

if nargin < 2
    mu = 0;
end
if nargin < 3
    beta = 1;
end

u = exp(-exp(-(x-mu)/beta));


end