function z = evaluate_rho(x, y, p, tau, q, cell_constraints)


max_n_cp = size(p,2) - 1;

z = cell(length(cell_constraints),max_n_cp+1);

rho = @(z,tau)(z.*(tau-(z<0)));

for i_constraint = 1:length(cell_constraints)
    for j_noc = 1:max_n_cp+1
        fv = zeros(length(tau),1);
        for k_quantile = 1:length(tau)
            if isempty(q{i_constraint,j_noc})
                curr_comp = 1; curr_quantile_index = k_quantile;
            else
                curr_comp = sum(q{i_constraint,j_noc} <= tau(k_quantile)) + 1;
                if curr_comp == 1
                    curr_quantile_index = k_quantile;
                else
                    curr_quantile_index = sum(tau(1:k_quantile) > q{i_constraint,j_noc}(curr_comp - 1));
                end
            end
            
            if i_constraint ~= length(cell_constraints)
                curr_alpha = p{i_constraint, j_noc}{curr_comp}(1);
                curr_beta = p{i_constraint, j_noc}{curr_comp}(2);
                curr_z = p{i_constraint, j_noc}{curr_comp}(2 + curr_quantile_index);
                
                fv(k_quantile) = sum(rho(y - x * curr_alpha - curr_z * x .^ curr_beta, tau(k_quantile)));
            else
                curr_alpha = p{i_constraint, j_noc}{curr_comp}(1);
                curr_c = p{i_constraint, j_noc}{curr_comp}(2);
                curr_beta = p{i_constraint, j_noc}{curr_comp}(3);
                curr_z = p{i_constraint, j_noc}{curr_comp}(3 + curr_quantile_index);
                
                fv(k_quantile) = sum(rho(y - curr_c - x * curr_alpha - curr_z * x .^ curr_beta, tau(k_quantile)));
            end
        end
        
        z{i_constraint, j_noc} = fv;
    end
end

% vec1 = n*log(tau.*(1-tau))' - n - n*log(Z{2,1}/n);
% vec2 = n*log(tau.*(1-tau))' - n - n*log(Z{2,2}/n);
% vec3 = n*log(tau.*(1-tau))' - n - n*log(Z{2,3}/n);
% vec4 = n*log(tau.*(1-tau))' - n - n*log(Z{2,4}/n);
% vec5 = n*log(tau.*(1-tau))' - n - n*log(Z{2,5}/n);
% Mat = [vec1, vec2, vec3, vec4, vec5];
% diff(Mat, [], 2)
% diff([n*log(tau.*(1-tau))' - n - n*log(Z{1,1}/n),n*log(tau.*(1-tau))' - n - n*log(Z{1,2}/n),n*log(tau.*(1-tau))' - n - n*log(Z{1,3}/n),n*log(tau.*(1-tau))' - n - n*log(Z{1,4}/n),n*log(tau.*(1-tau))' - n - n*log(Z{1,5}/n)],[],2)


end