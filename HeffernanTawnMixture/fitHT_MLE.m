function [Prm, NLOGL] = fitHT_MLE(HT, noc, plotOn, printOn)
% fit the Heffernan-Tawn mixture model using the MLE estimator -> found by an iterative algorithm

if nargin < 3
	plotOn = false;
end
if nargin < 4
	printOn = false;
end

% Initialise
eps = 1e-5;


alpha0 = -sort(-rand(1,noc)*0.9);
beta0 = zeros(1,noc);
mu0 = zeros(1,noc);
sig0 = ones(1,noc);

Rsd = zeros(length(HT.x),noc);
L = zeros(length(HT.x),noc);


alphaHat = alpha0;
betaHat = beta0;
muHat = mu0;
sigHat = sig0;
curr_cp = zeros(1,noc);

tx = HT.x; logtx = log(tx);
ty = HT.y;


cnt = 1; curr_diff = 10;

fun = @(p, LL)(HT_mixture_likelihood(p, tx, ty, LL, logtx));
while (cnt == 1 || curr_diff > eps) && cnt < 100
    if printOn
        fprintf('%.2f, %d\n',curr_diff,cnt)
    end
    tx_beta = exp(log(tx) .* betaHat(cnt,:));
    for i = 1:noc
        Rsd(:,i) = (ty - alphaHat(cnt,i) * tx - muHat(cnt,i) * tx_beta(:,i)) ./ (sigHat(cnt,i) * tx_beta(:,i));
        L(:,i) = normpdf(Rsd(:,i));
    end
    
    [newprm, newfval] = fminsearch(@(p)(fun(p, L)), [alphaHat(cnt,:); betaHat(cnt,:); muHat(cnt,:); sigHat(cnt,:)], optimset('MaxIter',3200*noc,'MaxFunEvals',4800*noc));
    
    alphaHat(cnt+1,:) = newprm(1,:);
    betaHat(cnt+1,:) = newprm(2,:);
    muHat(cnt+1,:) = newprm(3,:);
    sigHat(cnt+1,:) = newprm(4,:);
    curr_cp(cnt+1,:) = mean(L ./ sum(L,2) );
    
    curr_diff = sum(abs(alphaHat(cnt+1,:) - alphaHat(cnt,:)) ...
        + abs(betaHat(cnt+1,:) - betaHat(cnt,:)) ...
        + abs(muHat(cnt+1,:) - muHat(cnt,:)) ...
        + abs(sigHat(cnt+1,:) - sigHat(cnt,:)));
    
    cnt = cnt + 1;
end

if plotOn
	figure(1); clf; 
	subplot(2,2,1); plot(alphaHat); ylabel('\alpha');
	subplot(2,2,2); plot(betaHat); ylabel('\beta');
	subplot(2,2,3); plot(muHat); ylabel('\mu');
	subplot(2,2,4); plot(sigHat); ylabel('\sigma');
end

Prm = {newprm,curr_cp(end,:)};
NLOGL = newfval;

end

