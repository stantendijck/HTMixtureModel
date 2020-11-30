function NLOGL = GPD_like(x, mu, sigma, xi)

if xi < 0 && max(x) > mu - sigma/xi
    NLOGL = inf;
    return
end
if sigma < 0
    NLOGL = inf;
    return;
end

if min(x) < mu
    NLOGL = inf;
    return;
end

n = length(x);
NLOGL = n*log(sigma) + sum((1/xi + 1) * log(1+xi*(x-mu)/sigma));



end