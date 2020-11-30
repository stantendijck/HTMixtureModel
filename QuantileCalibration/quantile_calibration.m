function [output_q, output_p] = quantile_calibration(p, q, lower_middle_quantile, upper_middle_quantile, PltOn)

if nargin < 5
    PltOn = true;
end

%% Get the probabilities
p_lower = p(p <= lower_middle_quantile);
p_middle = p(p > lower_middle_quantile & p < upper_middle_quantile);
p_upper = p(p >= upper_middle_quantile);

%% Get the probabilities
q_lower = q(p <= lower_middle_quantile);
q_middle = q(p > lower_middle_quantile & p < upper_middle_quantile);
q_upper = q(p >= upper_middle_quantile);

%% Define where we want to calculate the function
mL = min(p_lower); ML = max(p_lower);
mM = min(p_middle); MM = max(p_middle);
mU = min(p_upper); MU = max(p_upper);

output_p_lower  = [exp(-(10:-1:0))*mL, mL:.005:ML];
output_p_middle = mM:.005:MM;
output_p_upper = [mU:.005:MU, 1 - exp(-(0:10)) * (1-MU)];


%% Approximate the empirical cdf
if length(p_lower) == 1
    debug = 1;
end
output_q_lower = quantile_cal_extrapolation(q_lower, p_lower, output_p_lower);
% output_p_middle = quantile_cal_interpolation(q_middle, p_middle, output_q_middle);
[output_q_middle, output_p_middle] = quantile_cal_np_interpolation(q_middle, p_middle, output_p_middle);
% output_p_middle = output_p_middle'; output_q_middle = output_q_middle';
output_q_upper = quantile_cal_extrapolation(q_upper, p_upper, output_p_upper);

if isempty(output_q_middle)
    output_q = [output_q_lower, output_q_upper];
    output_p = [output_p_lower, output_p_upper];
    return
end

%% Combine
if isempty(output_q_lower)
    I1 = [];
else
    I1 = find(output_q_middle < output_q_lower(end));
end
if isempty(output_q_upper)
    I2 = [];
else
    I2 = find(output_q_middle > output_q_upper(1),1);
end

indices_q = [1, length(output_q_middle)];
if ~isempty(I1)
    indices_q(1) = I1(end) + 1;
end
if ~isempty(I2)
    indices_q(2) = I2 - 1;
end

if isempty(output_p_lower)
    I1 = [];
else
    I1 = find(output_p_middle < output_p_lower(end));
end
if isempty(output_p_upper)
    I2 = [];
else
    I2 = find(output_p_middle > output_p_upper(1),1);
end



indices_p = [1, length(output_p_middle)];
if ~isempty(I1)
    indices_p(1) = I1(end) + 1;
end
if ~isempty(I2)
    indices_p(2) = I2 - 1;
end

indices = [max(indices_p(1),indices_q(1)),min(indices_p(2),indices_q(2))];

if any(indices == 0)
    I = find(output_p_upper > output_p_lower(end),1);
    output_q = [output_q_lower, output_q_upper(I:end)];
    output_p = [output_p_lower, output_p_upper(I:end)];
    return
end

%% Interpolate between (lower and middle) and (middle and upper) parts
N = 10;
boundsq = cell(2,2);
boundsq{1,1} = output_q_lower(end);
boundsq{1,2} = output_q_middle(indices(1));
boundsq{2,1} = output_q_middle(indices(2));
boundsq{2,2} = output_q_upper(1);

boundsq{1,1} = boundsq{1,1} + ((boundsq{1,2}-boundsq{1,1})/(N-1));
boundsq{1,2} = boundsq{1,2} - ((boundsq{1,2}-boundsq{1,1})/(N-1));
boundsq{2,1} = boundsq{2,1} + ((boundsq{2,2}-boundsq{2,1})/(N-1));
boundsq{2,2} = boundsq{2,2} - ((boundsq{2,2}-boundsq{2,1})/(N-1));

intra_q = cell(2,1);
intra_q{1} = boundsq{1,1}:((boundsq{1,2}-boundsq{1,1})/(N-1)):boundsq{1,2};
intra_q{2} = boundsq{2,1}:((boundsq{2,2}-boundsq{2,1})/(N-1)):boundsq{2,2};

boundsp = cell(2,2);
boundsp{1,1} = output_p_lower(end);
boundsp{1,2} = output_p_middle(indices(1));
boundsp{2,1} = output_p_middle(indices(2));
boundsp{2,2} = output_p_upper(1);

boundsp{1,1} = boundsp{1,1} + ((boundsp{1,2}-boundsp{1,1})/(N-1));
boundsp{1,2} = boundsp{1,2} - ((boundsp{1,2}-boundsp{1,1})/(N-1));
boundsp{2,1} = boundsp{2,1} + ((boundsp{2,2}-boundsp{2,1})/(N-1));
boundsp{2,2} = boundsp{2,2} - ((boundsp{2,2}-boundsp{2,1})/(N-1));

intra_p = cell(2,1);
intra_p{1} = boundsp{1,1}:((boundsp{1,2}-boundsp{1,1})/(N-1)):boundsp{1,2};
intra_p{2} = boundsp{2,1}:((boundsp{2,2}-boundsp{2,1})/(N-1)):boundsp{2,2};

output_q = [output_q_lower, intra_q{1}, output_q_middle(indices(1):indices(2)), intra_q{2}, output_q_upper];
output_p = [output_p_lower, intra_p{1}, output_p_middle(indices(1):indices(2)), intra_p{2}, output_p_upper];

%% Plot results
if PltOn
    figure; clf;
    plot(output_q_middle, output_p_middle,'LineWidth',4);
    hold on;
    if ~isempty(output_q_upper)
        plot(output_q_upper, output_p_upper,'LineWidth',4);
        plot(q_upper, p_upper, 'ko');
    end
    if ~isempty(output_q_upper)
        plot(output_q_lower, output_p_lower,'LineWidth',4);
        plot(q_lower, p_lower, 'ko');
    end
    plot(q_middle, p_middle,'k*');

    plot(output_q, output_p,'k--','LineWidth',2)
    xlim([output_q_lower(11) output_q_upper(end-11)])
end


