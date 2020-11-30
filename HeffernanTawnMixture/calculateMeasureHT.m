function pB = calculateMeasureHT(HT, xl, xu, yl, yu, rng_seed)
%% Calculate the probability of extreme set B for whole MCMC chain HT.p
% B = {(x,y): xl < x < xu; yl < y < yu}
rng(rng_seed);
rng_seed2 = randi(10000);
rng_seed3 = randi(10000);

%% Ge the randomness of the computation
rng(rng_seed2);
U_allocate_Rsd = rand(length(HT.x), 1);

% Model parameters
N = 1000; Delta = 1;
while xu - xl > 10
    xu = xu - 1;
end

if abs(yu) > abs(yl)
    while yu - yl > 10
        yu = yu - 1;
    end
else
    while yu - yl > 10
        yl = yl + 1;
    end
end
    
while Delta > 1/100
    N = 2*N;
    Delta = (xu - xl)/N;
end

%% Loop over MCMC chain
S = zeros(HT.max_iter, 1);
logHTx = log(HT.x);
for iter = 1:HT.max_iter
    %% Get parameters
    curr_prm = HT.p{iter};
    
    %% Get data allocated to component ic
    HTx_beta = exp(logHTx .* curr_prm(2,:));
    Rsd = getRsd(HT, curr_prm(1,:), HTx_beta, zeros(size(curr_prm(1,:))), ones(size(curr_prm(1,:))));
    for ic = 1:HT.noc
        Rsd(:,ic) = (Rsd(:,ic) - curr_prm(3,ic)) / curr_prm(4,ic);
    end
    L = getL(Rsd, HTx_beta, curr_prm(4,:));
    allocate = getallocate(L, U_allocate_Rsd);
    
    %% Loop over fits of the components
    for ic = 1:HT.noc
        %% Get parameters of component ic
        alpha = curr_prm(1, ic);
        beta = curr_prm(2, ic);
        
        % xic = HT.x(HT.allocate{iter, ic});
        % yic = HT.y(HT.allocate{iter, ic});
        
        xic = HT.x(allocate{ic});
        yic = HT.y(allocate{ic});
        
        %% Get residuals of component ic
        zic = (yic - alpha * xic) ./ (exp(log(xic) .* beta));
        
        %% Compute integral int_(xl)^(xu) int_(yl)^(yu) dP
        xvec = xl + ((1:N) - 1/2) * Delta;
        yvec = alpha * xvec + exp(log(xvec) .* beta) .* zic;
        S(iter) = S(iter) + Delta * sum(exp(-xvec)/2 .* sum(yl < yvec & yvec < yu, 1));
%         figure(101); hold on; plot(xvec,yvec,'k.');
    end
%     plot(HT.x,HT.y,'r.','MarkerSize',12);
end

pB = S / length(HT.x);

end