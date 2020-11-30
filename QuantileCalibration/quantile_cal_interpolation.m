function output_p = quantile_cal_interpolation(input_q, input_p, output_q)

% d_input_p = diff(input_p);
% sigma = max(d_input_p(d_input_p ~= 0)) * 2;
% kernel = @(x, y)(normcdf(x, y, sigma));
% 
% 
% output_p = zeros(size(output_q));
% for i_output_q = 1:length(output_q)
%     S = 0;
%     for i_input_q = 1:length(input_q)
%         S = S + kernel(output_q(i_output_q), input_q(i_input_q));
%     end
%     output_p(i_output_q) = S/length(input_q) * range(input_p) + input_p(1);
% end

nrminv = @(p,z)(norminv(z,p(1),p(2)));
calibration_fun = @(p)(sum((nrminv(p(1:2),input_p) - input_q).^2));
opt_prm = fminsearch(calibration_fun,[input_q(round(end/2)),1]);
% calibration_fun = @(p)(sum(((nrminv(p(1:2),input_p)+nrminv(p(3:4),input_p))/2 - input_q).^2));
% opt_prm = fminsearch(calibration_fun,[0,1,0,1]);
% 
% figure;
% plot(input_q, input_p, 'k*');
% hold on;
% plot(output_q, normcdf(output_q, opt_prm(1), opt_prm(2)));
% 
output_p = normcdf(output_q, opt_prm(1), opt_prm(2));
% output_p = (normcdf(output_q, opt_prm(1), opt_prm(2)) + normcdf(output_q, opt_prm(3), opt_prm(4)))/2;



% d_input_q = diff(input_q);
% d_input_p = diff(input_p);
% d_output_q = diff(output_q);
% 
% sigma = 2 * min(d_input_q);
% kernel = @(x, y)(normpdf(x, y, sigma));
% output_p = zeros(size(output_q));
% 
% for i_output_q = 1:length(d_output_q)
%     S = 0;
%     for i_input_q = 1:length(input_q)
%         S = S + kernel(output_q(i_output_q), input_q(i_input_q)) * input_p(i_input_p);
%     end
%     output_p(i_output_q) = S/length(input_q);
% end








end