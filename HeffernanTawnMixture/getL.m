function L = getL(Rsd, HTx_beta, sigma_curr)
%Get likelihood given residuals

L = zeros(size(Rsd));
for i = 1:size(Rsd,2)
    L(:,i) = 1/sqrt(2*pi) * 1./(sigma_curr(i) * HTx_beta(:,i)) ...
        .* exp( -Rsd(:,i).^2/2 );
end

end