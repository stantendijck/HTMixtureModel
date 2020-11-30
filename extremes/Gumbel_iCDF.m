function x = Gumbel_iCDF(u, mu, beta)
%Compute the inverse of the cumulative distribution function of a Gumbel
%random variable with parameters (mu, beta) -> default = (mu = 0, beta = 1)

if nargin < 2
    mu = 0;
end
if nargin < 3
    beta = 1;
end

x = mu - beta * log(-log(u));

end