function [TP_laplace, HS_laplace, GPD_prm, HS_original, TP_original] = ...
    PIT2(stormpeak, gpd_prm, orgHS, orgTP)

%% PIT transform the observations between xl and xu in Drc
u_Tp = zeros(length(stormpeak(:,3)),1);
u_HS = zeros(length(stormpeak(:,1)),1);

qTP = gpd_prm(1,4);
qHS = gpd_prm(2,4);

I_Tp = stormpeak(:,3) > quantile(stormpeak(:,1),qTP);
I_HS = stormpeak(:,1) > quantile(stormpeak(:,1),qHS);

prm_Tp = gpd_prm(1,:);
prm_HS = gpd_prm(2,:);

u_HS(~I_HS) = epit(stormpeak(~I_HS,1),orgHS(orgHS<quantile(orgHS,qHS)))*qHS;
u_HS(I_HS) = GPD_CDF(stormpeak(I_HS,1),prm_HS(1),prm_HS(2),prm_HS(3))*(1-qHS) + qHS;

u_Tp(~I_Tp) = epit(stormpeak(~I_Tp,3),orgTP(orgTP<quantile(orgTP,qTP)))*qTP;
u_Tp(I_Tp) = GPD_CDF(stormpeak(I_Tp,3),prm_Tp(1),prm_Tp(2),prm_Tp(3))*(1-qTP) + qTP;

Laplace_HS = Laplace_iCDF(u_HS);
Laplace_Tp = Laplace_iCDF(u_Tp);


TP_laplace = Laplace_Tp;
HS_laplace = Laplace_HS;
GPD_prm = [prm_Tp;prm_HS];
HS_original = stormpeak(:,1);
TP_original = stormpeak(:,3);

end
