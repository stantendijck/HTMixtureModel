function HT = fitHT_special(HT, alpha, beta, max_iter, keef)
% This fits the HT mixture model but it fixes alpha and beta and estimates
% optimal mu and sigma

%% Initialisation
rng(12345);

iter = 1;
n_cp = length(alpha) - 1;
noc = n_cp + 1;

burnIn = 1000;
AdpItr = 1000;

%% Pre-allocation
p = cell(max_iter, 1);                      % parameter estimates
new_p = zeros(2*noc,0);                     % parameter estimates in matrix
allocate = cell(max_iter, 1);               % allocation vector

%% Starting value
p{1} = [zeros(1,n_cp+1);ones(1,n_cp+1)];
new_p(:,1) = reshape(p{1},2*noc,1);

%% Step-size
stepsize.h = 0.05 * ones(2,n_cp+1);
stepsize.h1 = 0.05 * ones(2,n_cp+1);

%% MCMC
HT.logx = log(HT.x);
HT.accept_vec = zeros(max_iter,1);
while iter <= max_iter
    %% Print progress
    if mod(iter, max_iter) == 0; fprintf('\n');
    elseif mod(iter, 1000) == 0; fprintf('+');
    end
    
    %% Get current parameter values
    curr_prm = p{iter};
    
    mu_curr = curr_prm(1,:);
    sigma_curr = curr_prm(2,:);
    
    %% Calculate residuals and likelihoods
    HTx_beta = exp(HT.logx .* beta);
    Rsd = getRsd(HT, alpha, HTx_beta, mu_curr, sigma_curr);
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
    if iter <= burnIn + AdpItr % fixed proposal for the initial iterations
        proposed = propose(curr_prm, 1, stepsize);
    else % adaptive MCMC proposal
        proposed = propose(curr_prm, 2, stepsize);
    end
    
    %% Accept the new parameters
    t_curr_prm = [alpha;beta;curr_prm];
    t_proposed_prm = [alpha;beta;proposed];
    [t_curr_prm, HT] = accept(2, iter, HT, HTx_beta, t_curr_prm, allocate, t_proposed_prm, keef);
    curr_prm = t_curr_prm(3:4,:);
    
    %% Save the new (resp. old) parameters if accepted (rejected)
    p{iter+1} = curr_prm;
    new_p(:,iter+1) = reshape(p{iter+1}, 2*noc, 1);
    
    %% Next iteration
    iter = iter + 1;
end

HT.x = HT.x;
HT.y = HT.y;

HT.max_iter = max_iter;
HT.n_cp = n_cp;
HT.noc = noc;

HT.p = p;
