function fval = evaluate_beta(b, x, y, tau)
    
tol = 1e-4;
x_b = x.^b;

[alpha_begin, alpha_end] = find_feasible_alpha_grid(b, x, y, x_b);

a_hat = GoldenSectionSearch_ALPHA(@(a)(opt_rho_HT([a;beta],x, y, tau, 1, x_b)), alpha_begin, alpha_end, tol);
f = opt_rho_HT([a_hat;b],x, y, tau, 2, x_b);
fval = f{1};

end