function [output_q, output_p, cps] = Calibrate_Data(QR, j_noc, PltOn)

%% Get upper and lower quantile functions
[lower_quantile_funs, upper_quantile_funs] = get_upper_and_lower_quantile_functions(QR, j_noc);

%% Get some grid of quantiles from the model
[components, ~, ~] = get_quantiles(QR, lower_quantile_funs, upper_quantile_funs, j_noc);

%% Quantile calibration
cps = [0; QR.q{1, j_noc}; 1];

output_q = cell(j_noc, 1); output_p = cell(j_noc, 1);
for ic = 1:j_noc
    %     fprintf('Quantile Calibration on component %d/%d\n',ic, j_noc)
    rescaled_quantiles = (components{ic}{2} - cps(ic)) / (cps(ic + 1) - cps(ic));
    if 1
        [output_q{ic}, output_p{ic}] = quantile_calibration(rescaled_quantiles', components{ic}{1}', .1, .9, PltOn);
    else
        output_q{ic} = components{ic}{1}';
        output_p{ic} = rescaled_quantiles';
    end
end

if PltOn
    % Plot
    figure(7); clf;
    for ic = 1:j_noc
        plot(components{ic}{1}, components{ic}{2},'r*');
        hold on;
        plot(output_q{ic}, output_p{ic} * (cps(ic + 1) - cps(ic)) + cps(ic) ,'k-', 'LineWidth',2);
    end
    xlim([-10 10]);
end
