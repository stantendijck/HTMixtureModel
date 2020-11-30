function [alpha_begin, alpha_end] = find_feasible_alpha_grid(beta, x, y, x_beta)
% Define a grid alpha_grid = [a1, a2, a3] such that (ai, beta) are feasible
% parameter choices according to the Keef constraints

if nargin < 4
    x_beta = x.^beta;
end

alpha_mid = 0;
if ~all(HeffernanTawn.KeefConstraints(x,y,[alpha_mid;beta],[0;1],x_beta))
    alpha_mid_1 = alpha_mid - 0.1;
    while alpha_mid_1 > -.99 && ~all(HeffernanTawn.KeefConstraints(x,y,[alpha_mid_1 - 0.1;beta],[0;1],x_beta))
        alpha_mid_1 = alpha_mid_1 - 0.1;
    end
    alpha_mid_2 = alpha_mid + 0.1;
    while alpha_mid_2 < .99 && ~all(HeffernanTawn.KeefConstraints(x,y,[alpha_mid_2 + 0.1;beta],[0;1],x_beta))
        alpha_mid_2 = alpha_mid_2 + 0.1;
    end
    if alpha_mid_1 < - .99 && alpha_mid_2 > .99
        alpha_begin = [];
        alpha_end = [];
        return;
    end
    if alpha_mid_1 > - .99
        alpha_end = alpha_mid_1;
        alpha_begin = alpha_mid_1 - 1;
    elseif alpha_mid_2 < .99
        alpha_begin = alpha_mid_2;
        alpha_end = alpha_mid_2 + 0.1;
    end
else
    alpha_begin = alpha_mid;
    while alpha_begin > - 0.99 && all(HeffernanTawn.KeefConstraints(x,y,[alpha_begin - 0.1;beta],[0;1],x_beta))
        alpha_begin = alpha_begin - 0.1;
    end
    alpha_end = alpha_mid;
    while alpha_end < 0.99 && all(HeffernanTawn.KeefConstraints(x,y,[alpha_end + 0.1;beta],[0;1],x_beta))
        alpha_end = alpha_end + 0.1;
    end
end




