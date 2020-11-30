function vec_org = transform2original(vec_lap, GPD_prm, data_org)


quantiles_laplace = Laplace_CDF(vec_lap);
I_tp_below_thr = quantiles_laplace < GPD_prm(4);

vec_org = zeros(size(vec_lap));

if ~all(~I_tp_below_thr)
    vec_org(I_tp_below_thr) = quantile(data_org,quantiles_laplace(I_tp_below_thr));
end
if ~all(I_tp_below_thr)
    vec_org(~I_tp_below_thr) = GPD_iCDF((quantiles_laplace(~I_tp_below_thr) - GPD_prm(4))/(1 - GPD_prm(4)),GPD_prm(3),GPD_prm(1),GPD_prm(2));
end


end
