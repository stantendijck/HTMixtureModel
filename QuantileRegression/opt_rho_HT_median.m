function Y = opt_rho_HT_median(prm, x, y, tau, n_out, x_beta)
% prm = [alpha,beta], we penalise abs(alpha) < 1; beta <= 1 is implied by
% implementation method

if 1
    %Get the parameters
    alpha = prm(1);
    beta = prm(2);
    
    if length(tau) > 1
        begin_tau = max(0,3/2*tau(1) - tau(2)/2);
        end_tau = min(1,3/2*tau(end) - tau(end-1)/2);
        median_tau = (begin_tau + end_tau)/2;
    else
        median_tau = tau;
    end
    
    c_hat = GoldenSectionSearch(@(c)(sum((y-alpha*x - c) .* (median_tau - ((y - alpha*x - c) < 0)))), -10, 10, 1e-4);
%     tic()
%     c_hat = myfmin(@(c)(sum((y-alpha*x - c) .* (median_tau - ((y - alpha*x - c) < 0)))), 0);
%     toc()
%     tic()
%     c_hat = fminsearch(@(c)(sum((y-alpha*x - c) .* (median_tau - ((y - alpha*x - c) < 0)))), 0);
%     toc
    %calculate x^beta
    if nargin < 5
        x_beta = exp(log(x) .* beta);
    end
    
    %Check Keef constraints
    factor = 1 + 100*(~all(HeffernanTawn.KeefConstraints(x,y,[alpha;beta],[0,1],x_beta)));
    
    %initialise
    tol = 1e-4;
    opt_z = zeros(length(tau),1);
    fval = zeros(length(tau),1);
    z_start = -10; z_end = 10;
    tau_start = tau(1); tau_end = tau(end);
    
    %solve for tau(1)
    fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c_hat)) .* (tau_start - (y - (alpha*x + z*x_beta + c_hat)<0)))));
    opt_z(1) = GoldenSectionSearch_Z(fun, z_start, z_end, tol);
    
    %solve for tau(end)
    if length(tau) > 1
        fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c_hat)) .* (tau_end - (y - (alpha*x + z*x_beta + c_hat)<0)))));
        opt_z(end) = GoldenSectionSearch_Z(fun, opt_z(1), opt_z(1) + z_end - z_start, tol);
    end
    
    %solve for tau(i)
    if length(tau) > 2
        for i = 2:length(tau)-1
            fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c_hat)) .* (tau(i) - (y - (alpha*x + z*x_beta + c_hat)<0)))));
            opt_z(i) = GoldenSectionSearch_Z(fun, opt_z(i-1), opt_z(end), tol);
        end
    end
    
    %get the function values
    for i = 1:length(tau)
        fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c_hat)) .* (tau(i) - (y - (alpha*x + z*x_beta + c_hat)<0)))));
        fval(i) = fun(opt_z(i));
    end
    
    
    if length(tau) <= 2 || n_out == 2
        w = ones(1,length(tau));
    else
        w = diff([tau_start,(tau(2:end) + tau(1:end-1))/2,tau_end])/min(diff(tau));
    end
    
    if n_out == 1
        Y = factor * (w * fval);
    else
        Y = {factor * (w * fval), c_hat, opt_z};
    end
    
else
    % Get the parameters
    alpha = prm(1);
    beta = prm(2);
    c = prm(3);

    % calculate x^beta
    if nargin < 6
        x_beta = exp(log(x) .* beta);
    end

    % Check Keef constraints
    factor = 1 + 100*(~all(HeffernanTawn.KeefConstraints(x,y,[alpha;beta],[0,1],x_beta)));

    % initialise
    tol = 1e-4;
    opt_z = zeros(length(tau),1);
    fval = zeros(length(tau),1);
    z_start = -10; z_end = 10;
    tau_start = tau(1); tau_end = tau(end);

    % solve for tau(1)
    fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c)) .* (tau_start - (y - (alpha*x + z*x_beta + c)<0)))));
    opt_z(1) = GoldenSectionSearch_Z(fun, z_start, z_end, tol);

    % solve for tau(end)
    if length(tau) > 1
        fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c)) .* (tau_end - (y - (alpha*x + z*x_beta + c)<0)))));
        opt_z(end) = GoldenSectionSearch_Z(fun, opt_z(1), opt_z(1) + z_end - z_start, tol);
    end

    % solve for tau(i)
    if length(tau) > 2
        for i = 2:length(tau)-1
            fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c)) .* (tau(i) - (y - (alpha*x + z*x_beta + c)<0)))));
            opt_z(i) = GoldenSectionSearch_Z(fun, opt_z(i-1), opt_z(end), tol);
        end
    end

    % get the function values
    for i = 1:length(tau)
        fun = @(z)(sum(abs((y - (alpha*x + z*x_beta + c)) .* (tau(i) - (y - (alpha*x + z*x_beta + c)<0)))));
        fval(i) = fun(opt_z(i));
    end


    if length(tau) <= 2 || n_out == 2
        w = ones(1,length(tau));
    else
        w = diff([tau_start,(tau(2:end) + tau(1:end-1))/2,tau_end])/min(diff(tau));
        w = w/max(w);
    end

    if n_out == 1
        Y = factor * (w * fval);
    else
        Y = {factor * (w * fval), opt_z};
    end
end


end