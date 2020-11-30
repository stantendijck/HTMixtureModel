function x_min = GridSearch_BETA(fun, x0, x1, noPrint, tol)

beta_grid = [x0,x1];
while beta_grid(2)-beta_grid(1) > tol
    beta_grid = linspace(x0,x1,min(21,2+round((x1-x0)/tol)));
    
    fprintf('Grid size = %d, Grid diff = %.2f, tol = %.2f\n',length(beta_grid), beta_grid(2) - beta_grid(1), tol);
    fval = zeros(length(beta_grid),1);

    for i_beta = 1:length(beta_grid)
        fval(i_beta) = fun(beta_grid(i_beta));
    end
    fval

    [~,I] = min(fval);
    if I == length(beta_grid) || I == 1
        x_min = beta_grid(I);
        return
    end
    
    x0 = beta_grid(I-1);
    x1 = beta_grid(I+1);
end

[~,I] = min(fval);
x_min = beta_grid(I);
    



end