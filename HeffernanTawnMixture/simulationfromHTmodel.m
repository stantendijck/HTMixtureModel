function [X, Y] = simulationfromHTmodel(n,HT_K,HT,fiddle,TP_laplace,HS_laplace)
% input: cal: calsea, calswell, calall
% input: HT: HTsea, HTswell, HTall

X = zeros(n,1); Y = zeros(n,1);

% Simulate the indices of exceedances of the HT threshold
thr = Laplace_iCDF(fiddle.HT_quant);
U = rand(n,1);
Ithr = U < fiddle.HT_quant;
nbelow = sum(Ithr);
nexc = sum(~Ithr);

% Simulate from empirical model | X < thr
I = TP_laplace < thr;
population = [TP_laplace(I),HS_laplace(I)];
Ind = randsample(sum(I),nbelow,true);
X(Ithr) = population(Ind,1);
Y(Ithr) = population(Ind,2);



% First calculate mixture probabilities
% HT{HT_K}.p

curr_prm = HT{HT_K}.p{round(end/2)};

alpha_curr = curr_prm(1,:);
beta_curr = curr_prm(2,:);
mu_curr = curr_prm(3,:);
sigma_curr = curr_prm(4,:);

%% Calculate residuals and likelihoods
HTx_beta = exp(HT{HT_K}.logx .* beta_curr);
Rsd = getRsd(HT{HT_K}, alpha_curr, HTx_beta, mu_curr, sigma_curr);
L = getL(Rsd, HTx_beta, sigma_curr);
allocate = getallocate(L);
allocate_prob = zeros(HT_K,1);
for i = 1:HT_K
    allocate_prob(i) = mean(L(:,i) ./ sum(L,2));
end
allocate_prob_cumsum = [0;cumsum(allocate_prob)]';


% Simulate from QR model(1) | X > thr
thr = Laplace_iCDF(fiddle.HT_quant);
condX = Laplace_iCDF(U(~Ithr));

% Simulate components
component = sum(allocate_prob_cumsum < rand(nexc,1),2);

% allocate
Ysim = zeros(size(condX));

for i = 1:HT_K
    I = component == i;
    alp = alpha_curr(i);
    bet = beta_curr(i);
    m = mu_curr(i);
    s = sigma_curr(i);
    Ysim(I) = alp * condX(I) + condX(I) .^ bet .* (m + s * randsample(Rsd(allocate{i}==1,i),sum(I),true));
end

X(~Ithr) = condX;
Y(~Ithr) = Ysim;

% figure; plot(Xsim,Ysim,'k.');

