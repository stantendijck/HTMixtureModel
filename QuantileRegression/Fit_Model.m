function [p, f, t, q] = Fit_Model(x, y, tau, max_n_cp, model_constraint, noPrint)

if nargin < 5
    noPrint = false;
end

if max_n_cp >= length(tau)
    error('Fit_Model() does not work for max_n_cp > len(tau)');
end

if 1
    cell_constraints = {{'HT','median'}};
else
    cell_constraints = {{'HT'},{'HT','median'}};
end

temp_p = cell(length(cell_constraints), 1);
f = zeros(length(cell_constraints), max_n_cp + 1);
t = cell(length(cell_constraints), 1);

p = cell(length(cell_constraints), max_n_cp + 1);
q = cell(length(cell_constraints), max_n_cp + 1);

for i = 1:length(cell_constraints)
    tic();
    [temp_p{i}, f(i,:), t{i}] = quantreg_2(x, y, tau, max_n_cp, noPrint, cell_constraints{i}, model_constraint);
    toc();
    for j = 1:max_n_cp+1
        p{i,j} = temp_p{i}{j}{1}; if j == 1; p{i,j} = {p{i,j}}; end
        % q{i,j} = temp_p{i}{j}{2};
        q{i, j} = get_mixture_probability(x, y, p{i,j}, temp_p{i}{j}{2}, t{i}, cell_constraints{i});
    end
end

end