function HT = fitHT(HT, max_iter, n_cp, omega)


%% Initialisation
% rng(12345);

iter = 1;
noc = 1;

burnIn = 1000;
AdpItr = 1000;

%% Pre-allocation
p = cell(max_iter, 1);                      % parameter estimates
new_p = zeros(4*noc, 0);                     % parameter estimates in matrix
allocate = cell(max_iter, 1);				% allocation vector

%% Starting values
alpha0 = .5*ones(1,n_cp+1);
if n_cp ~= 0
  alpha0(2:end) = (.1 + .8/n_cp):(.8/n_cp):.9;
end
p{1} = [alpha0(end:-1:1);zeros(2,n_cp+1);ones(1,n_cp+1)];
new_p(:,1) = reshape(p{1},4*noc,1);

%% Step-size
stepsize.h = 0.05 * ones(4,n_cp+1);
stepsize.h1 = 0.05 * ones(4,n_cp+1);

% Prior on alpha is uniform on [-1,1] with mass omega at 1.
% Density at [-1,1) = 1/2*1/(1-omega));
% Prob mass at 1 = omega;
prior = @(x)(1/2 * (1-omega) * (x >= -1 & x < 1) + omega * (x==1));


%% MCMC
HT.logx = log(HT.x);
HT.accept_vec = zeros(max_iter,1);
HT.Like = zeros(max_iter,1);
while iter <= max_iter
    %% Print progress
    if mod(iter, max_iter) == 0; fprintf('\n');
    elseif mod(iter, 1000) == 0; fprintf('+');
    end
    
    %% Get current parameter values
	curr_prm = p{iter};
	
    alpha_curr = curr_prm(1,:);
    beta_curr = curr_prm(2,:);
    mu_curr = curr_prm(3,:);
    sigma_curr = curr_prm(4,:);
    
    %% Calculate residuals and likelihoods
    HTx_beta = exp(HT.logx .* beta_curr);
    Rsd = getRsd(HT, alpha_curr, HTx_beta, mu_curr, sigma_curr);
    L = getL(Rsd, HTx_beta, sigma_curr);
    allocate{iter} = getallocate(L);
    
	
	%% Calculate/Update covariance matrix
    if iter == burnIn + AdpItr + 1 % we don't need it before anyway
        stepsize.mu_upd = mean(new_p(:,burnIn:end),2);
        stepsize.adp_h = cov(new_p(:,burnIn:end)');
    elseif iter > burnIn + AdpItr + 1
        if mod(iter,1000) == 0 % just making sure, the update does not suffer from numerical errors
            stepsize.mu_upd = mean(new_p(:,burnIn:end),2);
            stepsize.adp_h = cov(new_p(:,burnIn:end)');
        else
            stepsize.mu_upd = stepsize.mu_upd + 1/(iter-burnIn+1) * (new_p(:,iter) - stepsize.mu_upd);
            stepsize.adp_h = (iter-burnIn)/(iter-burnIn+1)*stepsize.adp_h + (iter-burnIn)/(iter-burnIn+1)^2*(new_p(:,iter) - stepsize.mu_upd)*(new_p(:,iter)-stepsize.mu_upd)';
        end
    end
	
	
    %% Propose new parameters
    curr_prm = p{iter};
    if iter <= AdpItr + burnIn
        proposed = propose(curr_prm, 1, stepsize);
        alphastep = stepsize.h(1);
    else
        proposed = propose(curr_prm, 2, stepsize);
        alphastep = stepsize.adp_h(1,1);
    end
    proposed(1) = min(1, proposed(1));  % Truncate alpha at 1
    
    cdf = @(z)(0.5 * erfc(-z ./ sqrt(2)));
    from_a_to_b = @(a, b)((1 - cdf((1-a)/alphastep))*(b==1) + 1/alphastep * normpdf((b-a)/alphastep)*(b~=1));
    
    % The update is centering a Gaussian around previous value with mass at 1
    
    
    %% Check if the HT conditions |alpha| <= 1, beta < 1, sigma >= 0 hold
    if abs(proposed(1)) > 1 || proposed(2) > 1 || proposed(4) < 0
        flag_prm = false;
    else
        flag_prm = true;
    end

    %% Accept or Reject proposed parameters
	[curr_prm, HT] = accept(3, iter, HT, HTx_beta, curr_prm, allocate, proposed, {prior,from_a_to_b,omega});
    
    %% Save the updated parameters
    p{iter+1} = curr_prm;
    new_p(:,iter+1) = reshape(p{iter+1}, 4*noc, 1);
    
    %% Next iteration
    iter = iter + 1;
end

HT.x = HT.x;
HT.y = HT.y;

HT.max_iter = max_iter;
HT.n_cp = n_cp;
HT.noc = noc;

HT.p = p;
