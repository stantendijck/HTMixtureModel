function vec_laplace = transform2laplace(vec_org, GPD_prm, data_org)

vec_laplace = zeros(length(vec_org),1);
for k = 1:length(vec_org)
    if vec_org(k) > GPD_prm(1,1)  % if the observation is larger than mu_GPD
        u = GPD_CDF(vec_org(k),GPD_prm(1,1),GPD_prm(1,2),GPD_prm(1,3));
        u = u * (1 - GPD_prm(1,4)) + GPD_prm(1,4);
        vec_laplace(k) = Laplace_iCDF(u);
    else
        vec_laplace(k) = Laplace_iCDF(mean(data_org < vec_org(k)));
    end
end

end
