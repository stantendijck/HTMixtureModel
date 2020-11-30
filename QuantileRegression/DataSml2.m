function [x, y] = DataSml2(model, n, p, PltOn)

if nargin < 4
    PltOn = false;
end

t_p = 0.9;
new_n = round(n/(1-t_p));
if strcmp(model,'log')
    [A_laplace, ~, ~] = simulate_alg21(new_n,2,0.5);
    x = A_laplace(:, 1);
    y = A_laplace(:, 2);
    % figure(1); clf; plot(x,y,'k.');
elseif strcmp(model,'alog')
    [A_laplace, ~, ~] = simulate_alg22_quicker(new_n, 2, 0.5, [p, 0.5, 1 - p, 0.5]);
    x = A_laplace(:, 1);
    y = A_laplace(:, 2);
    % figure(1); clf; plot(x,y,'k.');
elseif strcmp(model,'normal')
    alpha_norm = 0.75;
    
    Z = mvnrnd([0,0],[1,alpha_norm;alpha_norm,1],new_n);
    x = Laplace_iCDF(normcdf(Z(:,1)));
    y = Laplace_iCDF(normcdf(Z(:,2)));
elseif strcmp(model,'normal_log_mixture')
    [A_laplace, ~, ~] = simulate_alg21(new_n, 2, 0.5);
    % [A_laplace, ~, ~] = simulate_alg22_quicker(n, 2, 0.5, [0.5, 0.5, 0.5, 0.5]);
    u = A_laplace(:,1);
    v = A_laplace(:,2);
    
    laplace_cdf = @(x)(1/2 + 1/2 * sign(x) .* (1 - exp(-abs(x))));
    laplace_icdf = @(u)( - sign(u - 0.5) .* log(1 - 2*abs(u - 0.5)) );
    frechet_cdf = @(x)( exp(-(1./x)) );
    frechet_icdf = @(u)( (-log(u)).^(-1) );
    
    frechet_u = frechet_icdf(laplace_cdf(u));
    frechet_v = frechet_icdf(laplace_cdf(v));
    
    alpha_norm = 0.5;
    Z = mvnrnd([0,0],[1,alpha_norm;alpha_norm,1],new_n);
    frechet_Z = frechet_icdf(normcdf(Z));
    
    frechet_X = max(p*frechet_Z(:,1),(1-p)*frechet_u);
    frechet_Y = max(p*frechet_Z(:,2),(1-p)*frechet_v);
    
    x = laplace_icdf(frechet_cdf(frechet_X));
    y = laplace_icdf(frechet_cdf(frechet_Y));
else
    error('model name is not supported');
end

if PltOn
    figure(1); clf;
    plot(x,y,'k.');
end

end