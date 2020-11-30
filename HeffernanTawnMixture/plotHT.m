function plotHT(HT)

%% PLOTS
burnIn = 500;
p_mat = cell2mat(HT.p');

new_pmat = zeros(HT.noc*4,HT.max_iter+1);
for j = 1:4
    for ic = 1:HT.noc
        new_pmat(ic+HT.noc*(j-1),:) = p_mat(j,ic:HT.noc:end);
    end
end
new_pmat = new_pmat(:,burnIn + 1:end);

names = {'\alpha','\beta','\mu','\sigma'};
figure(1); clf;
for j = 1:4
    subplot(4,1,j); hold on;
    for ic = 1:HT.noc
        plot(new_pmat(ic+HT.noc*(j-1),:));
    end
	ylabel(names{j});
end
% savePics('figures/HT_mixture_fit_MCMC1.pdf',1,'poster');

new_names = cell(1, 4*HT.noc);
names1 = {'\alpha','\beta','\mu','\sigma'};
for i = 1:HT.noc
    for j = 1:4
        new_names{(j-1)*HT.noc+i} = sprintf('%s_%d',names1{j},i);
    end
end


if HT.noc == 1
    figure(2); clf; hold on;
    for k = 1:4
        for l = k:4
            subplot(4,4,(k-1)*4 + (l-1) + 1);
            plot(new_pmat(k,burnIn:end),new_pmat(l,burnIn:end),'k.');
            xlabel(new_names{k}); ylabel(new_names{l});
        end
    end
elseif HT.noc == 2
    figure(2); clf; hold on;
    for k = 1:8
        for l = k:8
            subplot(8,8,(k-1)*8 + (l-1) + 1);
            plot(new_pmat(k,burnIn:end),new_pmat(l,burnIn:end),'k.');
            xlabel(new_names{k}); ylabel(new_names{l});
        end
    end
end

figure(3); clf;
plot(HT.accept_rate);
xlabel('Iteration');
ylabel('Acceptance rate');
% Allocation_probabilities = cell(1,HT.noc);
% if HT.noc > 1
%     for ic = 1:HT.noc
%         Alloc_probabilities = zeros(length(HT.x),HT.max_iter);
%         for i = 1:HT.max_iter
%             Alloc_probabilities(:,i) = HT.L(:,ic,i)./sum(HT.L(:,:,i),2);
%         end
%         Allocation_probabilities{ic} = mean(Alloc_probabilities,2);
%     end
%     if HT.noc == 2
%         Allocation_probabilities{3} = Allocation_probabilities{2};
%         Allocation_probabilities{2} = zeros(length(HT.x),1);
%     end
% elseif HT.noc == 1
%     Allocation_probabilities{1} = ones(length(HT.x),1);
% end
% Allocation_probabilities = cell2mat(Allocation_probabilities);
% if HT.noc < 4
%     figure(7); clf;
%     scatter(HT.x,HT.y,5,Allocation_probabilities,'filled'); colorbar();
%     colormap jet;
%     title('Allocation probabilities');
%     % savePics('figures/HT_mixture_fit_MCMC3.pdf',1,'poster');
% end




