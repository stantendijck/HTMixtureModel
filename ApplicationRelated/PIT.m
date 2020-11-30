function [TP_laplace, HS_laplace, GPD_prm, HS_original, TP_original] = ...
    PIT(stormpeak, marg_quant)

xl = 182; xu = 265;

q_GPD = 0.505:0.005:0.98;
[~,i_qgpd] = min(abs(q_GPD-marg_quant));


%% PIT transform the observations between xl and xu in Drc
I = stormpeak(:,2) > xl & stormpeak(:,2) < xu;
data = stormpeak(I,:);

u_HS = zeros(length(data(:,1)),1);
u_Tp = zeros(length(data(:,3)),1);

q = q_GPD(i_qgpd);
GPD_threshold_HS = quantile(data(:,1),q);
GPD_threshold_Tp = quantile(data(:,3),q);

I_HS = data(:,1) > GPD_threshold_HS;
I_Tp = data(:,3) > GPD_threshold_Tp;

phat_HS = fminsearch(@(p)(GPD_like(data(I_HS,1), GPD_threshold_HS, p(1), p(2))), [1,0.1], optimset('Display','off'));
phat_Tp = fminsearch(@(p)(GPD_like(data(I_Tp,3), GPD_threshold_Tp, p(1), p(2))), [1,0.1], optimset('Display','off'));

prm_HS = [GPD_threshold_HS, phat_HS(1), phat_HS(2), q];
prm_Tp = [GPD_threshold_Tp, phat_Tp(1), phat_Tp(2), q];

u_HS(~I_HS) = epit(data(~I_HS,1)) * q;
u_HS(I_HS) = GPD_CDF(data(I_HS,1), prm_HS(1), prm_HS(2), prm_HS(3)) * (1-q) + q;

u_Tp(~I_Tp) = epit(data(~I_Tp,3)) * q;
u_Tp(I_Tp) = GPD_CDF(data(I_Tp,3), prm_Tp(1), prm_Tp(2), prm_Tp(3)) * (1-q) + q;

HS_laplace = Laplace_iCDF(u_HS);
TP_laplace = Laplace_iCDF(u_Tp);

GPD_prm = [prm_Tp; prm_HS];
HS_original = data(:,1);
TP_original = data(:,3);


end

