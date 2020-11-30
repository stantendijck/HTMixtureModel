function x = Frechet_iCDF(u, a, mu, sigma)
%Compute the inverse of the cumulative distribution function of a Frecehet
%random variable with parameters (a, mu, sigma) -> default = (a = 1, mu = 0, sigma = 1)

if nargin < 2
    a = 1;
end
if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

x = mu + sigma * (-log(u)) .^ (-1/a);

end