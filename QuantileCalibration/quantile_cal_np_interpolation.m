function [output_q, new_output_p] = quantile_cal_np_interpolation(input_q, input_p, output_p)
%Non parametric monotonic interpolation of quantiles, minimizing the sum of
%squares

if size(input_q,1) == 1
    input_q = input_q';
end
if size(input_p,1) == 1
    input_p = input_p';
end

[~,I] = sort(input_q);

begin_val = min(input_q);
end_val = max(input_q);

nKnots = 5;
eps = 1e-4;
dKnots = (end_val-begin_val+2*eps)/(nKnots-1);
knots = (begin_val-eps):dKnots:(end_val+eps);
slm = slmengine(input_q, input_p,'knots',knots,'increasing','on','Extrapolation','linear');
% plotslm(slm);

new_output_p = [];
output_q = [];
for i = 1:length(output_p)
    temp = slmeval(output_p(i),slm,-1);
    if ~isnan(temp)
        new_output_p(end+1) = output_p(i);
        output_q(end+1) = temp;
    end
end


% figure(2); clf;
% plot(input_q,input_p,'k*');
% hold on;
% plot(output_q, new_output_p,'c-','LineWidth',3)

end