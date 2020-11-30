function fval = evaluate_beta_HT_median(b, x, y, tau, optimising_rho, n_out)
    
tol = 1e-3;
x_b = x.^b;

[alpha_begin, alpha_end] = find_feasible_alpha_grid(b, x, y, x_b);

a_hat = GoldenSectionSearch_ALPHA(@(a)(optimising_rho([a;b],x, y, tau, 1, x_b)), alpha_begin, alpha_end, tol);
f = optimising_rho([a_hat;b],x, y, tau, 2, x_b);
if n_out == 1
    fval = f{1};
else
    fval = {f{1}, a_hat, f{2}, f{3}};
end

end