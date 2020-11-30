function [x, y] = simulate_alg11(n, d, alpha)
% simulate from symmetric logistic model using Algorithm 1.1 in Stephenson
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

if d > 3
    warning('In simulate_alg1, the simulation of Gamma(k) for k >= 4, is not the most efficient');
end

% pre-allocation
y = zeros(n,d);

p = get_p(d,alpha);
u = - log(rand(n,d));
T = u./sum(u,2);
for i = 1:n
    if d ~= 1
        k = randsample(d,1,true,p{d});
    else
        k = 1;
    end
    z = sum(-log(rand(k,1)));
    y(i,:) = 1./(z*T(i,:).^alpha);
end

x = log(y);
% figure(1);
% clf;
% plot(y(:,1),y(:,2),'k.');




end