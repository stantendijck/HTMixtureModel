function F = GPD_CDF(x, mu, sigma, xi)

x(x == mu) = mu + 1e-8;
F = 1 - (1 + xi * (x - mu)/sigma) .^(-1/xi);
I = x < mu;
F(I) = 0;
if xi < 0
    I = x > mu - sigma/xi;
    F(I) = 1;
end
    

end