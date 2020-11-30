function F = Frechet_CDF(x, a, mu, sigma)
%Compute the cumulative distribution function of a Frechet random variable
%with parameters (a, mu, sigma) -> default = (a = 1, mu = 0, sigma = 1)

if nargin < 2
    a = 1;
end
if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

F = exp(-((x-mu)/sigma).^(-a));

end