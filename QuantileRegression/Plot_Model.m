%% Plot

switch default.purpose
    case {'poster','beamer'}
        fontsize = 18;
        legend_flag = false;
    case 'paper'
        fontsize = 11;
        legend_flag = true;
end

% Legend for plot
legend_str = cell(length(model.tau)+1,1);
legend_str{1} = 'data';

cell_constraints = {'HT','HT + median'};

FigCnt = 2;
for i_constraint = 1:length(cell_constraints)
    for j_noc = 1:size(q,2)
        figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on; plot(x, y, 'k.')
        
        tx = linspace(min(x),max(x),1000);
        
        colors = jet(length(tau));
        
        for k_quantile = 1:length(t{i_constraint})
            if isempty(q{i_constraint,j_noc})
                curr_comp = 1; curr_quantile_index = k_quantile;
            else
                curr_comp = sum(q{i_constraint,j_noc} <= tau(k_quantile)) + 1;
                if curr_comp == 1
                    curr_quantile_index = k_quantile;
                else
                    curr_quantile_index = sum(tau(1:k_quantile) > q{i_constraint,j_noc}(curr_comp - 1));
                end
            end
            
            if i_constraint ~= length(cell_constraints)
                curr_alpha = p{i_constraint, j_noc}{curr_comp}(1);
                curr_beta = p{i_constraint, j_noc}{curr_comp}(2);
                curr_z = p{i_constraint, j_noc}{curr_comp}(2 + curr_quantile_index);
                
                plot(tx, curr_alpha * tx + curr_z * tx .^ curr_beta, 'LineWidth', 3, 'Color', colors(k_quantile,:))
                legend_str{k_quantile+1} = sprintf('%.2f*x + %.2f*x^{%.2f}', curr_alpha, curr_z, curr_beta);
                
            else
                curr_alpha = p{i_constraint, j_noc}{curr_comp}(1);
                curr_c = p{i_constraint, j_noc}{curr_comp}(2);
                curr_beta = p{i_constraint, j_noc}{curr_comp}(3);
                curr_z = p{i_constraint, j_noc}{curr_comp}(3 + curr_quantile_index);
                
                plot(tx, curr_c + curr_alpha * tx + curr_z * tx .^ curr_beta, 'LineWidth', 3, 'Color', colors(k_quantile,:))
                legend_str{k_quantile+1} = sprintf('%.2f + %.2f*x + %.2f*x^{%.2f}', curr_c, curr_alpha, curr_z, curr_beta);
            
            end
        end
        if legend_flag
            legend(legend_str,'Location','best','fontsize',fontsize/2);
        end
        xlabel('X','fontsize',fontsize);
        ylabel('Y','fontsize',fontsize);
        title(sprintf('Best fit with %d mixture(s)',j_noc),'fontsize',fontsize);
        if i_constraint == 2
            set(gca,'fontsize',fontsize);
            % savePics(sprintf('Figures/%s_%d_%s.pdf',model,j_noc,cell_constraints{i_constraint}),1,purpose)
            savePics(sprintf('figures/%s_fit%d.pdf',model.name,j_noc),default.saveOn,default.purpose)
        end
    end
end


% 
% 
% 
% 
% 
% if isfield(p,'opt_prm')
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on; plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau)
%         plot(tx, p.alpha(1) * tx + p.z(j) * tx .^ p.beta(1), 'LineWidth', 2, 'Color', colors(j,:))
%         legend_str{j+1} = sprintf('%.2f*x + %.2f*x^{%.2f}',p.alpha(1),p.z(j),p.beta(1));
%     end
%     legend(legend_str,'Location','best','fontsize',fontsize/2);
%     title('H0 using fminsearch + GS search (no HT constraints)');
% end
% 
% if isfield(p,'opt_prm_HT')
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on; plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau)
%         plot(tx, p.alpha_HT(1) * tx + p.z_HT(j) * tx .^ p.beta_HT(1), 'LineWidth', 2, 'Color', colors(j,:))
%         legend_str{j+1} = sprintf('%.2f*x + %.2f*x^{%.2f}',p.alpha_HT(1),p.z_HT(j),p.beta_HT(1));
%     end
%     legend(legend_str,'Location','best','fontsize',fontsize/2);
%     if saveOn
%         switch model
%             case 'log'
%                 model_nme = 'Logistic';
%             case 'alog'
%                 model_nme = 'Asymmetric logistic';
%             otherwise
%                 model_nme = model;
%         end
%         title(sprintf('%s model',model_nme),'fontsize',fontsize);
%         xlabel('X','fontsize',fontsize);
%         ylabel('Y','fontsize',fontsize);
%     else
%         title('H0 using fminsearch + GS search (with HT constraints)','fontsize',fontsize);
%     end
%     savePics(sprintf('figures/%s_with_no_cp',model), saveOn, purpose);
% end
% 
% if isfield(p,'opt_prm_HT_median')
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on; plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau)
%         plot(tx, p.c_HT_median(1) + p.alpha_HT_median(1) * tx + p.z_HT_median(j) * tx .^ p.beta_HT_median(1), 'LineWidth', 2, 'Color', colors(j,:))
%         legend_str{j+1} = sprintf('%.2f + %.2f*x + %.2f*x^{%.2f}',p.c_HT_median(1),p.alpha_HT_median(1) ,p.z_HT_median(j),p.beta_HT_median(1));
%     end
%     legend(legend_str,'Location','best','fontsize',fontsize/2);
%     if saveOn
%         switch model
%             case 'log'
%                 model_nme = 'Logistic';
%             case 'alog'
%                 model_nme = 'Asymmetric logistic';
%             otherwise
%                 model_nme = model;
%         end
%         title(sprintf('%s model',model_nme),'fontsize',fontsize);
%         xlabel('X','fontsize',fontsize);
%         ylabel('Y','fontsize',fontsize);
%     else
%         title('H0 using fminsearch + GS search (with HT constraints + median -> 0)','fontsize',fontsize);
%     end
%     savePics(sprintf('figures/%s_median_with_no_cp',model), saveOn, purpose);
% end
% 
% if isfield(p,'opt_prm_M')
%     colors = jet(length(q.tau_M));
% 
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on;
%     plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau_M)
%         cnt = sum(q.tau_M(j)>t.opt_cp(1));
%         plot(tx, p.alpha_M(cnt + 1)*tx + p.z_M(j)*tx.^p.beta_M(cnt + 1), 'LineWidth', 2, 'Color', colors(j,:));
%         legend_str{j+1} = sprintf('%.2f*x + %.2f*x^{%.2f}',p.alpha_M(cnt + 1),p.z_M(j),p.beta_M(cnt + 1));
%     end
%     legend(legend_str,'location','best','fontsize',fontsize/2);
%     I = find((q.tau_M(2:end) + q.tau_M(1:end-1))/2 == t.opt_cp(1));
%     title(sprintf('H1 using fminsearch + GS search with %.2f < tau cp < %.2f (no HT constraints)',q.tau_M(I),q.tau_M(I+1)),'fontsize',fontsize);
% end
% 
% if isfield(p,'opt_prm_HT_M')
%     colors = jet(length(q.tau_HT_M));
% 
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on;
%     plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau_HT_M)
%         cnt = sum(q.tau_HT_M(j)>t.opt_cp_HT(1));
%         plot(tx, p.alpha_HT_M(cnt + 1)*tx + p.z_HT_M(j)*tx.^p.beta_HT_M(cnt + 1), 'LineWidth', 2, 'Color', colors(j,:));
%         legend_str{j+1} = sprintf('%.2f*x + %.2f*x^{%.2f}',p.alpha_HT_M(cnt + 1),p.z_HT_M(j),p.beta_HT_M(cnt + 1));
%     end
%     legend(legend_str,'location','best','fontsize',fontsize/2);
%     I = find((q.tau_HT_M(2:end) + q.tau_HT_M(1:end-1))/2 == t.opt_cp_HT(1));
%     title(sprintf('H1 using fminsearch + GS search with %.2f < tau cp < %.2f (with HT constraints)',q.tau_HT_M(I),q.tau_HT_M(I+1)),'fontsize',fontsize);
% %     if saveOn
% %         switch model
% %             case 'log'
% %                 model_nme = 'Logistic';
% %             case 'alog'
% %                 model_nme = 'Asymmetric logistic';
% %             otherwise
% %                 model_nme = model;
% %         end
% %         title(sprintf('%s model',model_nme));
% %         xlabel('X');
% %         ylabel('Y');
% %     else
% %         title(sprintf('H1 using fminsearch + GS search with %.2f < tau cp < %.2f (with HT constraints)',q.tau_HT_M(I),q.tau_HT_M(I+1)));
% %     end
% %     savePics(sprintf('figures/%s_with_cp',model), saveOn, purpose);
% end
% 
% 
% 
% if isfield(p,'opt_prm_HT_M')
%     colors = jet(length(q.tau));
% 
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on;
%     plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau)
%         I = find(q.tau_HT_M == q.tau(j));
%         cnt = sum(q.tau(j)>t.opt_cp_HT(1));
%         plot(tx, p.alpha_HT_M(cnt + 1)*tx + p.z_HT_M(I)*tx.^p.beta_HT_M(cnt + 1), 'LineWidth', 2, 'Color', colors(j,:));
%         legend_str{j+1} = sprintf('%.2f*x + %.2f*x^{%.2f}',p.alpha_HT_M(cnt + 1),p.z_HT_M(I),p.beta_HT_M(cnt + 1));
%     end
%     legend(legend_str,'location','best','fontsize',fontsize/2);
%     
%     if saveOn
%         switch model
%             case 'log'
%                 model_nme = 'Logistic';
%             case 'alog'
%                 model_nme = 'Asymmetric logistic';
%             otherwise
%                 model_nme = model;
%         end
%         title(sprintf('%s model',model_nme),'fontsize',fontsize);
%         xlabel('X','fontsize',fontsize);
%         ylabel('Y','fontsize',fontsize);
%     else
%         I = find((q.tau_HT_M(2:end) + q.tau_HT_M(1:end-1))/2 == t.opt_cp_HT(1));
%         title(sprintf('H1 using fminsearch + GS search with %.2f < tau cp < %.2f (with HT constraints)',q.tau_HT_M(I),q.tau_HT_M(I+1)),'fontsize',fontsize);
%     end
%     savePics(sprintf('figures/%s_with_cp',model), saveOn, purpose);
% end
% 
% 
% 
% if isfield(p,'opt_prm_HT_median_M')
%     colors = jet(length(q.tau_HT_median_M));
% 
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on;
%     plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau_HT_median_M)
%         cnt = sum(q.tau_HT_median_M(j)>t.opt_cp_HT_median(1));
%         plot(tx, p.c_HT_median_M(cnt + 1) + p.alpha_HT_median_M(cnt + 1)*tx + p.z_HT_median_M(j)*tx.^p.beta_HT_median_M(cnt + 1), 'LineWidth', 2, 'Color', colors(j,:));
%         legend_str{j+1} = sprintf('%.2f + %.2f*x + %.2f*x^{%.2f}',p.c_HT_median_M(cnt + 1),p.alpha_HT_median_M(cnt + 1),p.z_HT_median_M(j),p.beta_HT_median_M(cnt + 1));
%     end
%     legend(legend_str,'location','best','fontsize',fontsize/2);
%     I = find((q.tau_HT_median_M(2:end) + q.tau_HT_median_M(1:end-1))/2 == t.opt_cp_HT_median(1));
%     title(sprintf('H1 using fminsearch + GS search with %.2f < tau cp < %.2f (with HT constraints)',q.tau_HT_median_M(I),q.tau_HT_median_M(I+1)),'fontsize',fontsize);
% %     if saveOn
% %         switch model
% %             case 'log'
% %                 model_nme = 'Logistic';
% %             case 'alog'
% %                 model_nme = 'Asymmetric logistic';
% %             otherwise
% %                 model_nme = model;
% %         end
% %         title(sprintf('%s model',model_nme));
% %         xlabel('X');
% %         ylabel('Y');
% %     else
% %         title(sprintf('H1 using fminsearch + GS search with %.2f < tau cp < %.2f (with HT constraints)',q.tau_HT_median_M(I),q.tau_HT_median_M(I+1)));
% %     end
% %     savePics(sprintf('figures/%s_median_with_cp',model), saveOn, purpose);
% end
% 
% if isfield(p,'opt_prm_HT_median_M')
%     colors = jet(length(q.tau));
% 
%     figure(FigCnt); FigCnt = FigCnt + 1; clf; hold on;
%     plot(x, y, 'k.')
%     tx = linspace(min(x),max(x),1000);
%     for j = 1:length(q.tau)
%         I = find(q.tau_HT_median_M == q.tau(j));
%         cnt = sum(q.tau(j)>t.opt_cp_HT_median(1));
%         plot(tx, p.c_HT_median_M(cnt + 1) + p.alpha_HT_median_M(cnt + 1)*tx + p.z_HT_median_M(I)*tx.^p.beta_HT_median_M(cnt + 1), 'LineWidth', 2, 'Color', colors(j,:));
%         legend_str{j+1} = sprintf('%.2f + %.2f*x + %.2f*x^{%.2f}',p.c_HT_median_M(cnt + 1),p.alpha_HT_median_M(cnt + 1),p.z_HT_median_M(I),p.beta_HT_median_M(cnt + 1));
%     end
%     legend(legend_str,'location','best','fontsize',fontsize/2);
%     if saveOn
%         switch model
%             case 'log'
%                 model_nme = 'Logistic';
%             case 'alog'
%                 model_nme = 'Asymmetric logistic';
%             otherwise
%                 model_nme = model;
%         end
%         title(sprintf('%s model',model_nme),'fontsize',fontsize);
%         xlabel('X','fontsize',fontsize);
%         ylabel('Y','fontsize',fontsize);
%     else
%         I = find((q.tau_HT_median_M(2:end) + q.tau_HT_median_M(1:end-1))/2 == t.opt_cp_HT_median(1));
%         title(sprintf('H1 using fminsearch + GS search with %.2f < tau cp < %.2f (with HT constraints)',q.tau_HT_median_M(I),q.tau_HT_median_M(I+1)),'fontsize',fontsize);
%     end
%     savePics(sprintf('figures/%s_median_with_cp',model), saveOn, purpose);
% %     savePics(sprintf('C:\\Users\\tendijck\\Documents\\PhD\\Main Phase\\Posters\\2020 STORi Conference\\Figures\\%s_median_with_cp',model), saveOn, purpose);
% end
% 
% 
% % if max_K >= 5
% %     figure(5); clf; hold on;
% %     plot(x, y, 'k.')
% %     tx = linspace(min(x),max(x),1000);
% %     for j = 1:length(tau)
% %         cnt = sum(tau(j)>t.opt_cp_nc(1));
% %         plot(tx, p.opt_prm_M_nc(1, cnt + 1)*tx + p.opt_prm_M_nc(1, 2*n_cp+j+2)*tx.^p.opt_prm_M_nc(1, n_cp+2+cnt), 'LineWidth', 2, 'Color', colors(j,:));
% %         legend_str{j+1} = sprintf('%.2f*x + %.2f*x^{%.2f}',p.opt_prm_M_nc(1, cnt + 1),p.opt_prm_M_nc(1, 2*n_cp+j+2),p.opt_prm_M_nc(1, n_cp+2+cnt));
% %     end
% %     legend(legend_str,'location','best');
% %     I = find((tau(2:end) + tau(1:end-1))/2 == t.opt_cp_nc(1));
% %     title(sprintf('Model 2 using fminsearch + GS search with %.2f < tau cp < %.2f',tau(I),tau(I+1)));
% % %     savePics(sprintf('figures/qr_on_%s',model));
% % end
% 
% % Plot normal quantiles versus tau -> to check for multimodality in
% % non-mixture model
% % figure(4); subplot(1,2,1); plot(tau,normcdf(p(3:end))); title('PP-plot');
% % subplot(1,2,2); plot(norminv(tau),p(3:end)); title('QQ-plot');
