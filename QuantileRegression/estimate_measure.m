function est = estimate_measure(thr, x, y, p, q, tau, PltOn)

% I = x > 3.5;
% x = x(I);
% y = y(I);

% A = {(x,y): y > f1(x), x > f2(x)}
% For easyness, we define f1(x) = c1 = 5, f2(x) = c2 = 5;

if nargin < 6
    PltOn = true;
end

f1 = @(x)(thr(2).*ones(1,length(x)));
f2 = @(y)(thr(1).*ones(1,length(y)));

j_noc = 3;

prm = p{2, j_noc};
cp = q{2, j_noc};

quantile_funs = cell(length(tau),1);
prms = cell(length(tau),1);

for k_quantile = 1:length(tau)
    if isempty(cp)
        curr_comp = 1; curr_quantile_index = k_quantile;
    else
        curr_comp = sum(cp <= tau(k_quantile)) + 1;
        if curr_comp == 1
            curr_quantile_index = k_quantile;
        else
            curr_quantile_index = sum(tau(1:k_quantile) > cp(curr_comp - 1));
        end
    end

    
    curr_alpha = prm{curr_comp}(1);
    curr_c = prm{curr_comp}(2);
    curr_beta = prm{curr_comp}(3);
    curr_z = prm{curr_comp}(3 + curr_quantile_index);

    quantile_funs{k_quantile} = @(x)(curr_c + curr_alpha * x + curr_z * x .^ curr_beta);
    prms{k_quantile} = {curr_c, curr_alpha, curr_beta, curr_z, tau(k_quantile)};
end

% Need to find z such that the quantile function goes through (x_L, y_L)
point = fminsearch(@(X)(sqrt((f1(X(1))-X(2)).^2 + (f2(X(2))-X(1)).^2)),[0,0]);
xL = point(1); yL = point(2);

k_quantile = 1;
while k_quantile <= length(tau) && quantile_funs{k_quantile}(xL) < yL
    k_quantile = k_quantile + 1;
end

if k_quantile == 1
    error('Code is not made (yet) for k_quantile == 1');
end

if k_quantile == length(tau) + 1
    error('Code is not made (yet) for k_quantile == length(tau) + 1 -> need implentation of extrapolation of Z');
end

flag = false;
for i = 1:3
    if prms{k_quantile-1}{i} ~= prms{k_quantile}{i}
        flag = true;
    end
end
if flag
    error('Code is not made yet to deal with stuff near change points');
end

c = prms{k_quantile}{1};
alpha = prms{k_quantile}{2};
beta = prms{k_quantile}{3};

%% Binary search on [tau(k_quantile - 1), tau(k_quantile)]
lower_limit = tau(k_quantile - 1);
upper_limit = tau(k_quantile);

rho = @(z,tau)(z.*(tau-(z<0)));
eps = 1e-7;
while upper_limit - lower_limit > eps
    % Get new quantile
    new_q = (lower_limit + upper_limit) / 2;
    
    % Find new z given quantile
    fun = @(z)(sum(rho(y - c - alpha*x - z * x.^beta, new_q)));
    new_z = fminsearch(fun,0);
    
    % Get corresponding quantile regression function
    qr = @(x)(c + alpha*x + new_z * x.^beta);
    if qr(xL) > yL
        upper_limit = new_q;
    else
        lower_limit = new_q;
    end
end

%% Check whether the above works
if PltOn
    figure(1); clf; xv = 4:.01:6;
    plot(xv,quantile_funs{k_quantile}(xv)); hold on;
    plot(xv,quantile_funs{k_quantile-1}(xv));
    plot(xv,qr(xv),'r--')
    plot(xL,yL,'k*')
end

%% We have found:
corner_quantile = new_q;
corner_quantile_prm = {c, alpha, beta, new_z};

% Calculate lower bound on estimate:
% We need to use: X ~ Laplace(1); {(x,y): x > xL, y > q(x)} subset A
survival_cdf = 1/2 * sign(xL) * exp(-abs(xL));

lb_measure = survival_cdf * (1 - corner_quantile);

%% Approximate the other part of the area
% First find out what would be neglible
P = lb_measure/10;
upper_bound_x = -log(2-2*(1-P));

j_quantile = 1;
while j_quantile <= length(tau) && quantile_funs{j_quantile}(upper_bound_x) < f1(upper_bound_x)
    j_quantile = j_quantile + 1;
end

lower_limit = tau(j_quantile - 1);
upper_limit = corner_quantile;

if lower_limit < cp
    lower_limit = cp;
end

% lower_limit = 0.01;
% upper_limit = 0.99;
%%
N = 200;
grid = log((exp(lower_limit) : (exp(upper_limit) - exp(lower_limit)) / (N - 1) : exp(upper_limit)))';
% grid = (lower_limit : (upper_limit - lower_limit) / (N - 1) : upper_limit)';
dgrid = diff(grid);

x_grid = zeros(N, 1);
new_z = zeros(N, 1);
qr = cell(N, 1);

xbeta = x.^beta;
for i_quantile = 1:length(grid)
    % Find new z given quantile
    val = 0; % Choose like 4.5 to get more sensible results
    new_x = x(x > val);
    new_xbeta = xbeta(x > val);
    new_y = y(x > val);
    
%     fun = @(z)(sum(rho(y - c - alpha*x - z * xbeta, grid(i_quantile))));
    fun = @(z)(sum(rho(new_y - c - alpha*new_x - z * new_xbeta, grid(i_quantile))));
    new_z(i_quantile) = fminsearch(fun,0);
    
    % Get corresponding quantile regression function
%     qr{i_quantile} = @(x)(c + alpha*x + new_z(i_quantile) * x.^beta);
    qr{i_quantile} = @(x)(c + alpha*x + new_z(i_quantile) * x.^beta);
    
    % Solve qr(x) - f1(x) == 0
    x_grid(i_quantile) = fminsearch(@(x)(abs(qr{i_quantile}(x) - f1(x))), xL);
    old_z = new_z;
end

survival_cdf = 1/2 * sign(x_grid) .* exp(-abs(x_grid));

est_lower = sum(dgrid .* survival_cdf(1:end-1));
est_upper = sum(dgrid .* survival_cdf(2:end));

%% Check for how large the absolute error/relative error is by discretizing
est = (est_lower + est_upper)/2 + lb_measure;


if PltOn
    fprintf('\n');
    fprintf('Discretisation leads to lower bound of %f and upper bound of %f\n',est_lower,est_upper)
    fprintf('Absolute error = %f and relative error (wrt lb_measure) = %f\n',est_upper - est_lower, (est_upper - est_lower)/est)

    %% Plot the calculation method
    figure(2); clf;

    % Plot data
    plot(x,y,'k.'); hold on;
    xx = xlim; yy = ylim;

    % Plot region A
    xbox = [xL,xL:0.01:(xx(2)+1),xx(2)+1,f2((yy(2)+1):-.01:yL)];
    ybox = [yL,f1(xL:0.01:(xx(2)+1)),yy(2)+1,(yy(2)+1):-.01:yL];
    patch(xbox,ybox,'black','FaceAlpha',0.25)


    xlim([xL - 1 xx(2)]);
    ylim([yL - 2 yy(2)]);

    colors = jet(length(grid));
    xvec = min(x):0.01:(xx(2)+1);
    for i_quantile = 1:length(grid)-1
        plot(xvec,qr{i_quantile}(xvec),'Color',colors(i_quantile,:));
        plot(x_grid(i_quantile),qr{i_quantile}(x_grid(i_quantile)),'k*');

        d1 = (xx(2) + 1 - x_grid(i_quantile))/99;
        xbox = [x_grid(i_quantile):d1:(xx(2)+1),xx(2)+1:-d1:x_grid(i_quantile),x_grid(i_quantile)];
        ybox = [qr{i_quantile}(x_grid(i_quantile):d1:(xx(2)+1)),qr{i_quantile+1}(xx(2)+1:-d1:x_grid(i_quantile)),qr{i_quantile}(x_grid(i_quantile))];
        patch(xbox,ybox,colors(i_quantile,:),'FaceAlpha',0.25,'EdgeColor',colors(i_quantile,:));
    end

    plot(xvec,qr{end}(xvec),'Color',colors(end,:));
    plot(x_grid(end),qr{end}(x_grid(end)),'k*');

    d1 = (xx(2) + 1 - x_grid(i_quantile))/99;
    xbox = [x_grid(end):d1:(xx(2)+1),xx(2)+1,x_grid(end),x_grid(end)];
    ybox = [qr{end}(x_grid(end):d1:(xx(2)+1)),yy(2)+1,yy(2)+1,qr{end}(x_grid(end))];
    patch(xbox,ybox,colors(end,:),'FaceAlpha',0.5);
    title('Lower bound');

    %% Plot the calculation method
    figure(3); clf;

    % Plot data
    plot(x,y,'k.'); hold on;
    xx = xlim; yy = ylim;

    % Plot region A
    xbox = [xL,xL:0.01:(xx(2)+1),xx(2)+1,f2((yy(2)+1):-.01:yL)];
    ybox = [yL,f1(xL:0.01:(xx(2)+1)),yy(2)+1,(yy(2)+1):-.01:yL];
    patch(xbox,ybox,'black','FaceAlpha',0.25)


    xlim([xL - 1 xx(2)]);
    ylim([yL - 2 yy(2)]);

    colors = jet(length(grid));
    xvec = min(x):0.01:(xx(2)+1);
    for i_quantile = 1:length(grid)-1
        plot(xvec,qr{i_quantile}(xvec),'Color',colors(i_quantile,:));
        plot(x_grid(i_quantile),qr{i_quantile}(x_grid(i_quantile)),'k*');

        d1 = (xx(2) + 1 - x_grid(i_quantile+1))/99;
        xbox = [x_grid(i_quantile+1):d1:(xx(2)+1),xx(2)+1:-d1:x_grid(i_quantile+1),x_grid(i_quantile+1)];
        ybox = [qr{i_quantile}(x_grid(i_quantile+1):d1:(xx(2)+1)),qr{i_quantile+1}(xx(2)+1:-d1:x_grid(i_quantile+1)),qr{i_quantile}(x_grid(i_quantile+1))];
        patch(xbox,ybox,colors(i_quantile,:),'FaceAlpha',0.25,'EdgeColor',colors(i_quantile,:));
    end
    plot(xvec,qr{end}(xvec),'Color',colors(end,:));
    plot(x_grid(end),qr{end}(x_grid(end)),'k*');

    d1 = (xx(2) + 1 - x_grid(end))/99;
    xbox = [x_grid(end):d1:(xx(2)+1),xx(2)+1,x_grid(end),x_grid(end)];
    ybox = [qr{end}(x_grid(end):d1:(xx(2)+1)),yy(2)+1,yy(2)+1,qr{end}(x_grid(end))];
    patch(xbox,ybox,colors(end,:),'FaceAlpha',0.5);
    title('Upper bound');
end

end

















