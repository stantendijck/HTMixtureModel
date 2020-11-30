function pB = calculateMeasureHT_sml(HT, xl, xu, yl, yu, rng_seed)
%% Calculate the probability of extreme set B for whole MCMC chain HT.p
% B = {(x,y): xl < x < xu; yl < y < yu}
% Actually, code also runs for xl = constant vector, xu, yl, yu vectors
MCMC = 10000;

rng(rng_seed);
rng_seed2 = randi(10000);
rng_seed3 = randi(10000);

%% Ge the randomness of the computation
rng(rng_seed2);
U_allocate_Rsd = rand(length(HT.x), 1);

rng(rng_seed3);
U_simulate_X = rand(MCMC,1);
U_simulate_Z = randsample(length(HT.x),MCMC,true);

x = -log(1-U_simulate_X) + xl(1);
logx = log(x);

% Allocate y
y = zeros(size(x));

%% Loop over MCMC chain
logHTx = log(HT.x);

if HT.max_iter > 10000
    itervec = round(linspace(1,HT.max_iter,10000));
else
    itervec = 1:HT.max_iter;
end

S = zeros(length(itervec), length(xl));
cntr = 0;
for iter = itervec
    cntr = cntr + 1;
    %% Simulate from the HT mixture model
    % Get the parameters
    alpha = HT.p{iter}(1,:);
    beta = HT.p{iter}(2,:);
    mu = HT.p{iter}(3,:);
    sig = HT.p{iter}(4,:);
    
    % Get the residuals and allocations
    HTx_beta = exp(logHTx .* beta);
    Rsd = getRsd(HT, alpha, HTx_beta, mu, sig);
    L = getL(Rsd, HTx_beta, sig);
    allocate = getallocate(L, U_allocate_Rsd);
    
    % Loop over all HT mixture components
    for ic = 1:HT.noc
        inds = allocate{ic}(U_simulate_Z);
        z_comp_ic = Rsd(U_simulate_Z(inds), ic);
        x_comp_ic = x(inds);
        y_comp_ic = alpha(ic) * x_comp_ic + exp(logx(inds) .* beta(ic)) ...
            .* (mu(ic) + sig(ic) * z_comp_ic);
        y(inds) = y_comp_ic;
    end
    
    %% Calculate the number of points in the area R = [xl, xu] x [yl, yu]
    for iA = 1:length(xl)
        S(cntr,iA) = sum(xl(iA) < x & x < xu(iA) & yl(iA) < y & y < yu(iA));
    end
end

pB = (S * exp(-xl(1))/2) / MCMC;

end