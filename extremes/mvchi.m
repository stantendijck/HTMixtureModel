function [p, u] = mvchi(x)
%We need x to be on standard uniform margins, then this function 
%calculates chi(u) = P(X_(1:m-1) > u | X_m > u) for different values of u

[n, ~] = size(x);
max_u = 1 - 10/n;
% max_u = 0.8;
min_u = 0.5;
n_u = 1000;

new_u = zeros(1,1,n_u);
new_u(:,:) = linspace(min_u,max_u,n_u);

new_x = repmat(x,1,1,n_u);
new_p = mean(all(new_x>new_u,2))./mean(new_x(:,end,:)>new_u);

p = reshape(new_p,[1,n_u]);
u = reshape(new_u,[1,n_u]);

end