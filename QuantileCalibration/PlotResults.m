function PlotResults(curr_QR, probs, quantiles, j_noc, xl, xu, yl, yu, indices_used_quantiles, used_components, curr_quantile_funs, intersections)

%% Plot the lower bound
% colors = jet(length(indices_used_quantiles));
colors = jet(1e5);
cnt = 1;
for ic = used_components
    figure(j_noc+1+cnt); clf; %% Plot data
    
    plot(curr_QR.x,curr_QR.y,'k.'); hold on;
    xx = xlim; yy = ylim;
    
    cnt = cnt + 1;
end
xvec = (xl - 1):0.01:(xu+10);
for iquantile = 2:length(quantiles)-2
    if sum(iquantile == indices_used_quantiles) == 0
        continue
    end
    curr_comp = sum(curr_QR.q{2,j_noc} <= probs(iquantile)) + 1;

    figure(j_noc+1 + find(used_components == curr_comp)); hold on;
    plot(xvec, curr_quantile_funs{iquantile}(xvec) ,'Color',colors(round(probs(iquantile)*1e5),:));
    plot(intersections(iquantile,[1,3]),intersections(iquantile,[2,4]),'k*');
    
    xlimits = intersections(iquantile,[1,3]);
    xlimits_vec = linspace(xlimits(1), xlimits(2), 100);
    xbox = [xlimits_vec,xlimits_vec(end),xlimits_vec(end:-1:1),xlimits_vec(1)];
    ybox = [curr_quantile_funs{iquantile}(xlimits_vec),...
        curr_quantile_funs{iquantile+1}(xlimits_vec(end)),...
        curr_quantile_funs{iquantile+1}(xlimits_vec(end:-1:1)),...
        curr_quantile_funs{iquantile}(xlimits_vec(1))];
    patch(xbox,ybox,colors(round(probs(iquantile)*1e5),:),'FaceAlpha',0.75,'EdgeColor',colors(round(probs(iquantile)*1e5),:));
    cnt = cnt + 1;
end

cnt = 1;
for ic = used_components
    figure(j_noc+1+cnt); cnt = cnt + 1;
    
    %% Plot data
    plot(curr_QR.x,curr_QR.y,'k.'); hold on;

    %% Plot region A
    xbox = [xl,xu,xu,xl,xl];
    ybox = [yl,yl,yu,yu,yl];
    patch(xbox,ybox,'black','FaceAlpha',0.25)
    xlim([xl - 1 xx(2)]);
    ylim([yl - 2 yy(2)]);
    
    colorbar()
end

title('Lower and upper bound');
