function [best_p,best_fval,t] = quantreg_1(x, y, t, noPrint, constraints, model_constraint)
%Fit a quantile regression model with at most one changepoint

%% Check input arguments + Defaults
% if nargin == 1 && length(x) == 4 && iscell(x)
%     t = x{3};
%     y = x{2};
%     x = x{1};
% elseif nargin == 2 && size(x,2) == 2 && iscell(y)
%     if length(y) == 3
%         constraints = y{3};
%     else
%         constraints = {};
%     end
%     noPrint = y{2};
%     t = y{1};
%     y = x(:,2);
%     x = x(:,1);
% elseif nargin < 4
if nargin < 4
    noPrint = false;
end
if nargin < 5
end
if nargin < 6
    model_constraint = {};
end

%% Split t up into different sets
T = cell(length(t),length(t));
for i = 1:length(t)
    T{1,i} = t(1:i);
end
for i = 2:length(t)
    T{i,length(t)} = t(i:end);
end

%% Solve for all parts
p = cell(length(t), length(t));
fval = inf(length(t), length(t));
tau = cell(length(t), length(t));

if isempty(model_constraint)
    vec1 = 1:length(t);
    vec2 = 2:length(t);
else
    switch model_constraint{1}
        case 'size'
            vec1 = model_constraint{2}:length(t);
            vec2 = 2:(length(t)+1-model_constraint{2});
    end
end

for i = vec1
    [p{1, i}, fval(1, i), tau{1, i}] = quantreg_0(x, y, T{1, i}, noPrint, constraints);
end

for i = vec2
    [p{i, length(t)}, fval(i, length(t)), tau{i, length(t)}] = quantreg_0(x, y, T{i,length(t)}, noPrint, constraints);
end

fvals = zeros(length(t),1);
for i = 1:length(t)-1
    fvals(i) = fval(1,i) + fval(i+1,end);
end
fvals(end) = fval(1,length(t));


[~, I_fvals] = min(fvals);

if I_fvals ~= length(t)
    best_p = cell(2,1);
    
    % Get all parameters
    best_p{1} = cell(2,1);
    best_p{1}{1} = p{1, I_fvals};
    best_p{1}{2} = p{I_fvals + 1, length(t)};
    
    % Get the corresponding changepoint
    best_p{2} = (tau{1,I_fvals}(end) + tau{I_fvals + 1, length(t)}(1))/2;
    
    best_fval = fval(1, I_fvals) + fval(I_fvals + 1, length(t));
else
    best_p = cell(2,1);
    
    % Get all parameters
    best_p{1} = cell(1,1);
    best_p{1}{1} = p{1, length(t)};
    
    % Get the corresponding changepoint
    best_p{2} = [];
    
    best_fval = fvals(I_fvals);
end









end