function [pB_median_after_burnIn, estimates_mat, A_original, QR, HT, cal] = DataAnalysisWrapper(stormpeak, fiddle, model, default)

B = fiddle.B;
marg_quant = fiddle.marg_quant;
HT_quant = fiddle.HT_quant;

%% Bootstrap procedure
if B ~= 1
    I = randsample(size(stormpeak,1),size(stormpeak,1),true);
else
    I = 1:size(stormpeak,1);
end
temp_stormpeak = stormpeak(I,:);

%% Probability integral transform to standard Laplace margins
[TP_laplace, HS_laplace, GPD_prm, HS_original, TP_original] = PIT(temp_stormpeak, marg_quant);

%% Rectangles design
xgrid_original_scale = 12:.75:17;
ygrid_original_scale = 3.1:.75:22;

xgrid_marginal_scale = transform2laplace(xgrid_original_scale, GPD_prm(1,:), TP_original);
ygrid_marginal_scale = transform2laplace(ygrid_original_scale, GPD_prm(2,:), HS_original);

A_original = getRectangles(xgrid_original_scale,ygrid_original_scale);
A_marginal = getRectangles(xgrid_marginal_scale,ygrid_marginal_scale);

%% Modelling
x = TP_laplace; tx = x(x > quantile(x, HT_quant));
y = HS_laplace; ty = y(x > quantile(x, HT_quant));

%% Quantile Regression Model
% Fit the QR model
fprintf('Fitting QR Model\n');
QR = QuantileRegression(model, default, false);
[p, f, ~, q] = Fit_Model(tx, ty, model.tau, model.max_n_cp, model.constraint, default.noPrint);

QR.x = tx;
QR.y = ty;
QR.p = p;
QR.f = f;
QR.q = q;

% Calibrating the QR model (possibly using smoothing) and estimate the probabilities
fprintf('Calibrating QR Model\n');
cal.output_q = cell(model.max_noc, 1);
cal.output_p = cell(model.max_noc, 1);
cal.cps = cell(model.max_noc, 1);

for j_noc = 1:model.max_noc
    fprintf('%d: ',j_noc');
    [cal.output_q{j_noc}, cal.output_p{j_noc}, cal.cps{j_noc}] = Calibrate_Data(QR, j_noc, default.PltOn & j_noc==2);
end

estimates_mat = cell(size(A_marginal,1),1);

tic();
parfor iset = 1:size(A_marginal,1)
    
    curr_set = A_marginal(iset,:);
    
    xl = curr_set(1);
    xu = curr_set(2);
    yl = curr_set(3);
    yu = curr_set(4);
    
    fprintf('|\n\b');
    estimates = cell(1,4);
    for j_noc = 1:model.max_noc
        [S1, S2, probs, quantiles, indices_used_quantiles, used_components, curr_quantile_funs, intersections] = Calculate_Measure(QR, cal.output_q{j_noc}, cal.output_p{j_noc}, cal.cps{j_noc}, j_noc, xl, xu, yl, yu);
        estimates{j_noc} = (S1 + S2)/2;
    end
    estimates_mat{iset} = cell2mat(estimates);
end
fprintf('\n');
toc();

%% Heffernan-Tawn model

% Fit model
fprintf('Fitting HT mixture model\n');
HT = cell(model.max_noc, 1);
        
max_iter = 100000;
for j_noc = 1:model.max_noc
    HT{j_noc}.x = QR.x;
    HT{j_noc}.y = QR.y;
    HT{j_noc} = fitHT(HT{j_noc}, max_iter, j_noc - 1);
end

% pmat = cell2mat(HT{2}.p');
% 
% burnIn = 4000;
% prm = cell(8);
% for i = 1:4
%     prm{i} = pmat(i,(burnIn+1):2:end);
%     prm{i+4} = pmat(i,(burnIn+2):2:end);
% end
% 
% labels = {'\alpha_1','\beta_1','\mu_1','\sigma_1','\alpha_2','\beta_2','\mu_2','\sigma_2'};
% figure;
% cnt = 1;
% for i = 1:8
%     for j = 1:8
%         subplot(8,8,cnt); cnt = cnt + 1;
%         plot(prm{i},prm{j},'k.');
%         xlabel(labels{i});
%         ylabel(labels{j});
%     end
% end


% Calculate the probabilities using the HT mixture model
fprintf('Calculating measures using HT mixture model\n');
pB = cell(model.max_noc, size(A_marginal,1));
pB_median = zeros(model.max_noc, size(A_marginal,1));

pB_median_after_burnIn = cell(size(A_marginal,1),1);
for iset = 1:size(A_marginal,1)
    pB_median_after_burnIn{iset} = zeros(model.max_noc,1);
end

vec2 = cell(model.max_noc,1);

tic();
parfor j_noc = 1:model.max_noc
    vec = cell(length(1:length(ygrid_marginal_scale)-1:size(A_marginal,1)),1);
    for igroup = 1:length(ygrid_marginal_scale)-1:size(A_marginal,1)

        curr_set = A_marginal(igroup:igroup+length(ygrid_marginal_scale)-2,:);

        xl = curr_set(:,1);
        xu = curr_set(:,2);
        yl = curr_set(:,3);
        yu = curr_set(:,4);
        fprintf('\b|\n');

        vec{igroup} = calculateMeasureHT_sml(HT{j_noc}, xl, xu, yl, yu, 1000);
    end
    vec2{j_noc} = vec;
end


for j_noc = 1:model.max_noc
    for igroup = 1:length(ygrid_marginal_scale)-1:size(A_marginal,1)
        for iset = igroup:igroup+length(ygrid_marginal_scale)-2
            pB{j_noc, iset} = vec2{j_noc}{igroup}(:,iset-igroup+1);
            pB_median(j_noc, iset) = median(vec2{j_noc}{igroup}(:,iset-igroup+1));
            pB_median_after_burnIn{iset}(j_noc) = median(vec2{j_noc}{igroup}(1000:end,iset-igroup+1));
        end
    end
end



end