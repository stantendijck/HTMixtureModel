function [p,u] = chi_bar(x)
%We need x to be on standard uniform margins, then this function 
%calculates chi(u) = P(X_1 > u | X_2 > u) for different values of u

n = size(x,1);
max_u = 1 - 10/n;
min_u = 0.5;
u = linspace(min_u,max_u,1000);
p = zeros(size(u));

for i = 1:length(u)
%     p(i) = mean(x(x(:,1)>u(i),2) > u(i));
    p(i) = 2*log(1-u(i))./log(1-2*u(i) + mean(x(:,1)<u(i) & x(:,2)<u(i))) - 1;
    if isnan(p(i))
        p(i) = 0;
    end
end


end