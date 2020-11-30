function x = GPD_iCDF(u, xi, mu, sigma)
%Compute the inverse of the cumulative distribution function of a GPD
%random variable with parameters (xi, mu, sigma) -> default = (mu = 0, sigma = 1)

if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

if xi == 0
    x = mu - sigma * log(1-u);
else
    x = mu + sigma * ((1-u).^(-xi) - 1)/xi;
end

end