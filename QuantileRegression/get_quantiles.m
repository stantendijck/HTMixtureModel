function [components, curr_z, quantiles] = get_quantiles(QR, lower_quantile_funs, upper_quantile_funs, j_noc)

%% Fit the model
rho = @(z, t)(z .* (t - (z < 0)));

Begin = 0.001;
End = 0.999;
N = 1000;
eps = 1e-4;


cps = [Begin;QR.q{1, j_noc};End];

quantiles = [];
for i = 1:length(cps)-1
    quantiles = [quantiles, (cps(i) + eps):((cps(i+1) - cps(i) - 2*eps)/(N-1)):(cps(i+1) - eps)];
end
curr_z = zeros(length(quantiles), 1);



components = cell(j_noc, 1);
for ic = 1:j_noc
    components{ic}{1} = [];
    components{ic}{2} = [];
end

for iquantile = 1:length(quantiles)
    %% Get component
    curr_comp = sum(QR.q{1,j_noc} < quantiles(iquantile)) + 1;
    if curr_comp == 4
        debug = 1;
    end
    
    %% Get parameters of component
    curr_alpha = QR.p{1,j_noc}{curr_comp}(1);
    curr_c = QR.p{1,j_noc}{curr_comp}(2);
    curr_beta = QR.p{1,j_noc}{curr_comp}(3);
    curr_x_beta = QR.x .^ curr_beta;
    
    %% Fit the model once naively for quantiles(iquantile)
    fun = @(z, t)(sum(rho(QR.y - curr_alpha * QR.x - curr_c - z * curr_x_beta, t)));
    curr_z(iquantile) = fminsearch(@(z)(fun(z, quantiles(iquantile))), 0);
    curr_quantile = @(x)(curr_c + curr_alpha * x + curr_z(iquantile) * x .^ curr_beta);
    
    %% Find intersection with component below and above
    new_threshold_lower = min(QR.x); new_threshold_upper = min(QR.x);
    if curr_comp ~= 1
        new_threshold_lower = fminbnd(@(x)(abs(lower_quantile_funs{curr_comp - 1}(x) - curr_quantile(x))), min(QR.x), quantile(QR.x, .99));
    end
    if curr_comp ~= j_noc
        new_threshold_upper = fminbnd(@(x)(abs(upper_quantile_funs{curr_comp + 1}(x) - curr_quantile(x))), min(QR.x), quantile(QR.x, .99));
    end
    
    new_threshold = max(new_threshold_lower, new_threshold_upper);
    new_threshold = min(new_threshold, quantile(QR.x,.9));
    %% Fit model again but above new_threshold
    if cps(curr_comp+1)-cps(curr_comp) > 0.1 && abs(new_threshold - min(QR.x)) > 1e-4
        I = QR.x > new_threshold;
        yI = QR.y(I);
        xI = QR.x(I);
        xbetaI = curr_x_beta(I);
        fun = @(z, t)(sum(rho(yI - curr_alpha * xI - curr_c - z .* xbetaI, t)));
        curr_z(iquantile) = fminsearch(@(z)(fun(z, quantiles(iquantile))), 0);
    end
    
    components{curr_comp}{1} = [components{curr_comp}{1};curr_z(iquantile)];
    components{curr_comp}{2} = [components{curr_comp}{2};quantiles(iquantile)];
end
