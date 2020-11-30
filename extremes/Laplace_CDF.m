function u = Laplace_CDF(x, mu, b)
%Compute the cumulative distribution function of a Laplace random variable
%with parameters (mu, b) -> default = (mu = 0, b = 1)

if nargin < 2
    mu = 0;
end
if nargin < 3
    b = 1;
end

u = 1/2 + 1/2 * sign(x - mu) .* (1 - exp(-abs(x-mu)/b));



end