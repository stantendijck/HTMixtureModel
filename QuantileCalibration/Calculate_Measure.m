function [S1, S2, probs, quantiles, indices_used_quantiles, used_components, curr_quantile_funs, intersections] = Calculate_Measure(QR, output_q, output_p, cps, j_noc, xl, xu, yl, yu, PltOn)

if nargin < 10
    PltOn = false;
end

%% Get all probabilities and quantiles used
probs = [];
quantiles = [];
for ic = 1:j_noc
    probs = [probs, output_p{ic} * (cps(ic + 1) - cps(ic)) + cps(ic)];
    quantiles = [quantiles, output_q{ic}];
end

%% Loop over grid of quantiles
Laplace_CDF = @(x)(1 - 1/2 * exp(-x));

curr_quantile_funs = cell(length(probs),1);

S_upper = 0;
S_lower = 0;

indices_used_quantiles = [];
used_components = [];
intersections = zeros(length(probs),4);

comp_indices = zeros(length(probs),1);
for iprobs = 1:length(probs)
    %% Get current component
    curr_comp = sum(QR.q{1,j_noc} <= probs(iprobs)) + 1;
    comp_indices(iprobs) = curr_comp;
    
    if curr_comp == 4
        debug = 1;
    end
    
    %% Get parameters
    curr_alpha = QR.p{1,j_noc}{curr_comp}(1);
    curr_c = QR.p{1,j_noc}{curr_comp}(2);
    curr_beta = QR.p{1,j_noc}{curr_comp}(3);
    
    %% Find the intersection points if they exist of the quantile function and A
    curr_quantile_funs{iprobs} = @(x)(curr_c + curr_alpha * x + quantiles(iprobs) * exp(log(x) .* curr_beta));
    
    y_xl = curr_quantile_funs{iprobs}(xl);
    y_xu = curr_quantile_funs{iprobs}(xu);
    
    first_crossing = []; second_crossing = [];
    if yl < y_xl && y_xl < yu
        first_crossing = [xl, y_xl];
    end
    if yl < y_xu && y_xu < yu
        second_crossing = [xu, y_xu];
    end
    
    eps = 1e-2;
    fun = curr_quantile_funs{iprobs};
    [x_yl, fval_l] = fminbnd(@(x)(abs(fun(x) - yl)), xl - eps, xu + eps);
    [x_yu, fval_u] = fminbnd(@(x)(abs(fun(x) - yu)), xl - eps, xu + eps);
    
    
    if xl < x_yl && x_yl < xu % && abs(fval_u) < 1
        first_crossing = [x_yl, yl];
    elseif xl < x_yl && x_yl < xu && abs(fval_u) > 1
        debug = 1;
    end
    if xl < x_yu && x_yu < xu % && abs(fval_l) < 1
        second_crossing = [x_yu, yu];
    elseif xl < x_yu && x_yu < xu && abs(fval_l) > 1
        debug = 1;
    end
    
    if ~isempty(first_crossing) && ~isempty(second_crossing)
        if iprobs > 1
            S_lower(iprobs) = (probs(iprobs) - probs(iprobs-1)) * (Laplace_CDF(second_crossing(1))-Laplace_CDF(first_crossing(1)));
        end
        if iprobs < length(probs)
            S_upper(iprobs) = (probs(iprobs+1) - probs(iprobs)) * (Laplace_CDF(second_crossing(1))-Laplace_CDF(first_crossing(1)));
        end
        
        indices_used_quantiles(end + 1) = iprobs;
        used_components(end + 1) = curr_comp;
        intersections(iprobs,:) = [first_crossing, second_crossing];
    end
    
end

used_components = unique(used_components);

S1 = sum(S_upper);
S2 = sum(S_lower);


if PltOn
    fprintf('\n');
    fprintf('Discretisation leads to lower bound of %f and upper bound of %f\n',S2,S1)
    fprintf('Absolute error = %f and relative error (wrt lb_measure) = %f\n',S1 - S2, (S1 - S2)/S2)

    %% Plot the calculation method
    figure(2); clf;

    % Plot data
    plot(QR.x,QR.y,'k.'); hold on;
    xx = xlim; yy = ylim;

    % Plot region A
    xbox = [xl,xu,xu,xl,xl];
    ybox = [yl,yl,yu,yu,yl];
    patch(xbox,ybox,'black','FaceAlpha',0.25)


    xlim([xl - 1 xx(2)]);
    ylim([yl - 2 yy(2)]);

    colors = jet(length(indices_used_quantiles));
    xvec = min(QR.x):0.01:(xx(2)+1);
    for i_quantile = 1:length(indices_used_quantiles)
        Ind = indices_used_quantiles(i_quantile);
        plot(xvec,curr_quantile_funs{Ind}(xvec),'Color',colors(i_quantile,:));
        crossing = intersections(Ind,:);
        plot(crossing([1,3]),crossing([2,4]),'k*');
        
        if i_quantile ~= length(indices_used_quantiles)
            l = (crossing(3) - crossing(1))/99;
            xbox = [crossing(1):l:crossing(3),crossing(3),crossing(3):-l:crossing(1),crossing(1)];
            ybox = [curr_quantile_funs{Ind}(crossing(1):l:crossing(3)),curr_quantile_funs{Ind+1}([crossing(3),crossing(3):-l:crossing(1)]),curr_quantile_funs{Ind}(crossing(1))];
            patch(xbox,ybox,colors(i_quantile,:),'FaceAlpha',0.25,'EdgeColor',colors(i_quantile,:));
        end
    end
%     title('');


    %% Plot the calculation method
    figure(3); clf;

    % Plot data
    plot(QR.x,QR.y,'k.'); hold on;
    xx = xlim; yy = ylim;

    % Plot region A
    xbox = [xl,xu,xu,xl,xl];
    ybox = [yl,yl,yu,yu,yl];
    patch(xbox,ybox,'black','FaceAlpha',0.25)


    xlim([xl - 1 xx(2)]);
    ylim([yl - 2 yy(2)]);

    colors = jet(length(indices_used_quantiles));
    xvec = min(QR.x):0.01:(xx(2)+1);
    for i_quantile = 1:length(indices_used_quantiles)
        Ind = indices_used_quantiles(i_quantile);
        plot(xvec,curr_quantile_funs{Ind}(xvec),'Color',colors(i_quantile,:));
        crossingNew = intersections(Ind,:);
        plot(crossingNew([1,3]),crossingNew([2,4]),'k*');

        if i_quantile > 1
            crossingOld = intersections(Ind-1,:);
            l = (crossingOld(3) - crossingOld(1))/99;
            xbox = [crossingNew(1):l:crossingNew(3),crossingNew(3),crossingNew(3):-l:crossingNew(1),crossingNew(1)];
            ybox = [curr_quantile_funs{Ind-1}(crossingNew(1):l:crossingNew(3)),curr_quantile_funs{Ind}([crossingNew(3),crossingNew(3):-l:crossingNew(1)]),curr_quantile_funs{Ind-1}(crossingNew(1))];
            patch(xbox,ybox,colors(i_quantile,:),'FaceAlpha',0.25,'EdgeColor',colors(i_quantile,:));
        end
    end
    title('Estimate 2');
end








end
