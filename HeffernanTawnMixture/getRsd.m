function Rsd = getRsd(HT, alpha, HTx_beta, mu, sigma)
%Get Residuals of HT object with the parameters (alpha, beta, mu, sigma)

Rsd = zeros(length(HT.x),length(alpha));
for i = 1:length(alpha)
    Rsd(:,i) = (HT.y - alpha(i) * HT.x - mu(i) * HTx_beta(:,i))./(sigma(i)* HTx_beta(:,i));
end

end