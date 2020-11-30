function x = Laplace_iCDF(u, mu, b)
%Compute the inverse of the cumulative distribution function of a Laplace
%random variable with parameters (mu, b) -> default = (mu = 0, b = 1)

if nargin < 2
    mu = 0;
end
if nargin < 3
    b = 1;
end

x = mu - b * sign(u - 0.5) .* log(1 - 2*abs(u - 0.5));

end