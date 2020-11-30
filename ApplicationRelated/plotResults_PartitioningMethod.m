function plotResults_Application(x, y, K, estimates_mat, A, saveOn)

if nargin < 6
	saveOn = false;
end

% QR(2)
p = length(x.sea)/(length(x.sea) + length(x.swell));

plot(x.sea,y.sea,'k.'); hold on; plot(x.swell,y.swell,'k.'); 

XX = cellfun(@(x,y)(x(K)*p + y(K)*(1-p)),estimates_mat.sea,estimates_mat.swell);
% m = min(XX(XX>0));
% M = min(1,max(XX));
m = 1e-07;
M = 1;
a = (M/m)^(1/255);

% M = 0.0112;

epsilon = 0.05;

% Model 1: QR(1) -> estimates_mat{iset}(1);
% col = [linspace(1,0,256)',linspace(1,0,256)',linspace(1,0,256)'];
col = colormap(flipud(hot));
for iset = 1:size(A, 1)
    
    curr_set = A(iset, :);
    
    xl = curr_set(1);
    xu = curr_set(2);
    yl = curr_set(3);
    yu = curr_set(4);
    
    if XX(iset) < 1e-16
        continue
    end
    current_estimate = XX(iset);
    current_color = round(max(0,log(current_estimate/m)/log(a)))+1;
    patch([xl,xu,xu,xl,xl],[yl,yl,yu,yu,yl],col(current_color,:),'FaceAlpha',0.9,'EdgeAlpha',0);

%     current_color = max(min(1,log(1  + (exp(1)-1) * max(0,(XX(iset)  - m)/(M-m)))),0);
    
%     patch([xl,xu,xu,xl,xl],[yl,yl,yu,yu,yl],'red','FaceAlpha',current_color,'EdgeAlpha',0);
%     text(xl+epsilon,yl/2+yu/2,sprintf('%s',num2str(current_estimate,2)));
end
plot(x.sea,y.sea,'k.'); hold on; plot(x.swell,y.swell,'k.'); 
% title(sprintf('HT(%d): Estimates of sets',K));
% xlim([0 10.5])

% c.Ticks = linspace(m,M,11);
% c.Limits = [m,M];
c = colormap(col);
cbh = colorbar;
cbh.Ticks = [1e-5,1e-3,1e-1];
% caxis([m M]);
caxis([m 1]);
set(gca,'ColorScale','log');

% xlim([0 10.5])

savePics('OutputFig1.pdf',saveOn,'paper');

