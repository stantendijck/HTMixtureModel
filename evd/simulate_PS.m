function S = simulate_PS(n,alpha)
% Simulate from a positive stable distribution with parameter alpha
% 
% INPUT:
% n         number of samples
% alpha     (1 x d) vector of parameters, with d the dimension of the model
%
% OUTPUT
% S         (n x d) sample from PS(alpha)

% dimension
d = size(alpha,2);

% pre allocation
S = zeros(n,d);

% U = [Unif(0,pi), Exp(1)]
U = rand(n,2,d);
U(:,1,:) = U(:,1,:)*pi;
U(:,2,:) = -log(U(:,2,:));
for i_d = 1:d
    t_alpha = alpha(i_d);
    if t_alpha == 1
        S(:,i_d) = 1;
    else
        a1 = sin((1-t_alpha)*U(:,1,i_d));
        a2 = (a1./U(:,2,i_d)).^((1-t_alpha)/t_alpha);
        a3 = sin(t_alpha*U(:,1,i_d));
        a4 = sin(U(:,1,i_d)).^(1/t_alpha);
        S(:,i_d) = a2 .* a3 ./ a4;
    end
end



end