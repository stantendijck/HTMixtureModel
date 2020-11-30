function [x_laplace, x_gumbel, x_frechet] = simulate_alg21(n, d, alpha)
% simulate from symmetric logistic model using Algorithm 2.1 in Stephenson
% (2003)
%
% INPUT:
% n         number of samples
% d         dimension
% alpha     dependence parameter
%
% OUTPUT:
% x         (n x d) sample on Gumbel scale
% y         (n x d) sample on Frechet scale

% pre-allocation
x_frechet = zeros(n,d);

% S ~ PS(alpha)
S = simulate_PS(n,alpha);

% U = [d-Exp(1)]
U = -log(rand(n,d));

% set y = ((s/Exp(1))^(alpha),...,(s/Exp(1))^(alpha));
x_frechet = zeros(n,d);
for id = 1:d
    x_frechet(:,id) = (S./U(:,id)).^(alpha);
end

x_gumbel = log(x_frechet);
% figure(1);
% clf;
% plot(y(:,1),y(:,2),'k.');

x_unif = exp(-x_frechet.^(-1));
x_laplace = -sign(x_unif - 1/2) .* log(1 - 2*abs(x_unif-1/2));




end