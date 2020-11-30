function q = get_mixture_probability(x, y, prm, q0, tau, constraints)
% This function returns the "optimal" mixture probability given the
% estimates of the parameters in prm and the data [x, y].

% %% Get the Koenker rho function
flag_HT_median_constraints = false;
rho = @(z, t)(abs(z.*(t - (z<0))));
switch length(constraints)
    case 1
        if strcmp(constraints,'HT')
            optimising_rho = @(z, prm, xbeta, x, y, t)(sum(rho(y - prm * x - z * xbeta, t)));
        end
    case 2
        if all(strcmp(constraints,'HT') + strcmp(constraints,'median'))
            optimising_rho = @(z, prm, xbeta, x, y, t)(sum(rho(y - prm(2) - prm(1) * x - z * xbeta, t)));
            flag_HT_median_constraints = true;
        end
end

%%

q = nan(size(q0));

ncp = length(q0);

for icp = 1:ncp
    % Initialisation, i.e., get upper and lower bound
    itau = 1; flag = true;
    while flag && itau <= length(tau)
        if tau(itau) > q0(icp)
            flag = false;
        else
            itau = itau + 1;
        end
    end
    search_grid = [tau(itau - 1), tau(itau)];
    
    flag_not_equal = true;
    while flag_not_equal && search_grid(2) - search_grid(1) > 1e-8
        % Update current value
        curr_val = mean(search_grid);
        
        % Evaluate the functions
       if flag_HT_median_constraints
           alpha1 = prm{icp}(1); alpha2 = prm{icp+1}(1);
           c1 = prm{icp}(2); c2 = prm{icp+1}(2);
           beta1 = prm{icp}(3); beta2 = prm{icp+1}(3);
           [~, fval1] = fminsearch(@(z)(optimising_rho(z, [alpha1, c1], x.^beta1, x, y, curr_val)), prm{icp}(end));
           [~, fval2] = fminsearch(@(z)(optimising_rho(z, [alpha2, c2], x.^beta2, x, y, curr_val)), prm{icp+1}(1));
       else
           alpha1 = prm{icp}(1); alpha2 = prm{icp+1}(1);
           beta1 = prm{icp}(2); beta2 = prm{icp+1}(2);
           [~, fval1] = fminsearch(@(z)(optimising_rho(z, alpha1, x.^beta1, x, y, curr_val)), prm{icp}(end));
           [~, fval2] = fminsearch(@(z)(optimising_rho(z, alpha2, x.^beta2, x, y, curr_val)), prm{icp+1}(1));
       end
       
       % Binary search
       if fval1 < fval2
           search_grid = [curr_val, search_grid(2)];
       elseif fval2 < fval1
           search_grid = [search_grid(1), curr_val];
       else
           flag_not_equal = false;
       end
    end
    
    q(icp) = mean(search_grid);
end








end