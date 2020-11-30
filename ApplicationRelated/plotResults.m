function plotResults_QR(x, y, K, estimates_mat, A, saveOn)

if nargin < 6
	saveOn = false;
end

% QR(2)
plot(x,y,'k.'); hold on;


XX = cellfun(@(x)(x(K)),estimates_mat);
% m = min(XX(XX>0));
% M = min(1,max(XX));
m = 1e-7;
% m = 1e-4;
% M = 0.0112;
M = 1;
a = (M/m)^(1/255);
epsilon = 0.05;

% Model 1: QR(1) -> estimates_mat{iset}(1);
% col = autumn; col = col(end:-1:1,:);
% col = [linspace(1,0,256)',linspace(1,0,256)',linspace(1,0,256)'];
col = colormap(flipud(hot));
for iset = 1:size(A, 1)
    
    curr_set = A(iset, :);
    
    xl = curr_set(1);
    xu = curr_set(2);
    yl = curr_set(3);
    yu = curr_set(4);
    
    if estimates_mat{iset}(K) < 1e-16
        continue
    end
    
    current_estimate = estimates_mat{iset}(K);
    current_color = round(max(0,log(current_estimate/m)/log(a)))+1;
    patch([xl,xu,xu,xl,xl],[yl,yl,yu,yu,yl],col(current_color,:),'EdgeAlpha',0,'FaceAlpha',0.9);
    
%     current_color = max(min(1,log(1  + (exp(1)-1) * max(0,(estimates_mat{iset}(K)  - m)/(M-m)))),0);
%     patch([xl,xu,xu,xl,xl],[yl,yl,yu,yu,yl],'red','FaceAlpha',current_color,'EdgeAlpha',0);
%     text(xl+epsilon,yl/2+yu/2,sprintf('%s',num2str(current_estimate,2)));
end
plot(x,y,'k.');
% title(sprintf('QR(%d): Estimates of sets',K));
% xlim([0 10.5])

c = colormap(col);
% c.Limits = [m,M];

cbh = colorbar;
cbh.Ticks = [1e-5,1e-3,1e-1];
% caxis([m M]);
caxis([m 1]);
set(gca,'ColorScale','log');


savePics('OutputFig1.pdf',saveOn,'paper');

