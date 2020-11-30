function [best_p, best_fval, t] = quantreg_2(x, y, t, max_n_cp, noPrint, constraints, model_constraint)
%Fit a quantile regression model with at most max_n_cp changepoints ->
%works also optimally for max_n_cp = 0, max_n_cp = 1

%% Check input arguments + Defaults
% if nargin == 1 && length(x) == 4 && iscell(x)
%     max_n_cp = x{4};
%     t = x{3};
%     y = x{2};
%     x = x{1};
% elseif nargin == 2 && size(x,2) == 2 && iscell(y)
%     if length(y) == 4
%         constraints = y{4};
%     else
%         constraints = {};
%     end
%     noPrint = y{3};
%     max_n_cp = y{2};
%     t = y{1};
%     y = x(:,2);
%     x = x(:,1);
% elseif nargin < 5

if nargin < 5
    noPrint = false;
end

if nargin < 6
    constraints = {};
end

if nargin < 7
    model_constraint = {};
end

if isempty(model_constraint)
    flag_model_constraint = false;
else
    flag_model_constraint = true;
end

%% Split t up into different sets
T = cell(length(t),length(t));
for j = 1:length(t)
    for i = 1:length(t)
        T{i,j} = t(i:j);
    end
end

%% Ask quicker programs to solve the problem if max_n_cp = 0, 1

% Get all optimal parameters per value for K (number of mixture components)
best_p = cell(max_n_cp+1,1);
best_fval  = inf(max_n_cp + 1,1);


if max_n_cp == 0
    [best_p, best_fval, t] = quantreg_0(x, y, t, noPrint, constraints);
    return
elseif max_n_cp == 1
    [best_p0, best_fval0, ~] = quantreg_0(x, y, t, noPrint, constraints);
    [best_p1, best_fval1, t1] = quantreg_1(x, y, t, noPrint, constraints, model_constraint);
    
    best_p{1} = cell(2,1);
    best_p{1}{1} = best_p0;
    best_p{1}{2} = [];
    
    best_p{2} = cell(2,1);
    best_p{2}{1} = best_p1{1};
    best_p{2}{2} = best_p1{2};
    
    best_fval(1) = best_fval0;
    best_fval(2) = best_fval1;
    
    t = t1;
    return
end

%% Solve for all parts
p = cell(length(t), length(t));
fval = inf(length(t), length(t));
tau = cell(length(t), length(t));

for i = 1:length(t)
    parfor j = 1:length(t)
        if j < i
            continue
        end
        if flag_model_constraint && strcmp(model_constraint{1},'size') && j - i + 1 < model_constraint{2}
            continue
        end
        [p{i, j}, fval(i, j), tau{i, j}] = quantreg_0(x, y, T{i, j}, noPrint, constraints);
    end
end


% Get all optimal parameters per value for K (number of mixture components)
% best_p = cell(max_n_cp+1,1);
% best_fval  = inf(max_n_cp + 1,1);
opt_t = cell(max_n_cp+1,1);


% K = 0 has to go separately
best_fval(1) = fval(1,length(t));
best_p{1} = {p{1,length(t)};[]};
opt_t{1} = tau{1,length(t)};


for K = 1:max_n_cp
    best_p{K+1} = cell(2,1);    
    
    combinations = combnk(1:(length(t)-1),K);
    
    for ith_combination = 1:size(combinations, 1)
        curr_fval = 0;
        curr_prm = cell(K,1);
        curr_tau = cell(K,1);
        for jth_row = 1:(size(combinations, 2)+1)
            if jth_row == 1
                curr_fval = curr_fval + fval(1,combinations(ith_combination,jth_row));
                curr_prm{jth_row} = p{1,combinations(ith_combination,jth_row)};
                curr_tau{jth_row} = T{1,combinations(ith_combination,jth_row)};
            elseif jth_row == K+1
                curr_fval = curr_fval + fval(combinations(ith_combination,jth_row-1) + 1, length(t));
                curr_prm{jth_row} = p{combinations(ith_combination,jth_row-1) + 1, length(t)};
                curr_tau{jth_row} = T{combinations(ith_combination,jth_row-1) + 1, length(t)};
            else
                curr_fval = curr_fval + fval(combinations(ith_combination,jth_row-1) + 1, combinations(ith_combination,jth_row));
                curr_prm{jth_row} = p{combinations(ith_combination,jth_row-1) + 1, combinations(ith_combination,jth_row)};
                curr_tau{jth_row} = T{combinations(ith_combination,jth_row-1) + 1, combinations(ith_combination,jth_row)};
            end
        end
        if curr_fval < best_fval(K+1)
            best_fval(K+1)  = curr_fval;
            best_p{K+1}{1} = curr_prm;
            best_fval(K+1) = curr_fval;
            opt_t{K+1} = curr_tau;
            
            best_p{K+1}{2} = nan(length(curr_tau)-1, 1);
            for i = 1:length(curr_tau)-1
                best_p{K+1}{2}(i) = (curr_tau{i}(end) + curr_tau{i+1}(1))/2;
            end
        end
    end
end













end