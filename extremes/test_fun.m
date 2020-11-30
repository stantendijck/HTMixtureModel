% X = DistributionTransformers.COPPERWT;
% X = DistributionTransformers.TOTALWT;
X = DistributionTransformers.TERTRATKVA;
X(isnan(X)) = [];


K_CV = 10;

new_X = X(X > quantile(X,.9));
new_X = X(X > 2.5);

n = length(new_X);
I = randsample(n,n);

vec = round(linspace(1,n,K_CV+1));

mu = min(new_X);
F = zeros(n,1);
for iCV = 1:10
    test_X = new_X(I(vec(iCV):vec(iCV+1)));
    train_X = new_X([I(1:vec(iCV)-1);I(vec(iCV)+1:end)]);
    
    p0 = [1,1];
    opt_p = fminsearch(@(p)(GPD_like(train_X, mu, p(1), p(2))), p0);
    
    sigma = opt_p(1);
    xi = opt_p(2);
    F(I(vec(iCV):vec(iCV+1))) = GPD_CDF(test_X,mu,sigma,xi);
end

figure; plot(log(1-real(F))); hold on;
% figure; plot(real(F)); hold on;
plot([1,n],[log(1/n),log(1/n)],'r--');
title('COPPERWT')

vec = unique(DistributionTransformers.MAKER);

n = length(DistributionTransformers.MAKER);
makers = [];
for i = 1:n
    if DistributionTransformers.MANUFDTE(i) > 20000000
        makers = [makers;DistributionTransformers.MAKER(i)];
    end
end
makers = unique(makers);

X = DistributionTransformers.Plant_no;
n = length(X);
location = [];
for i = 1:n
    location = [location;str2num(X{i}(1:2))];
end
location = unique(location);


cnt = 4;
for j = 1:length(location)
    curr_location = location(j);
    figure(cnt); clf; cnt = cnt + 1;
    for i = 1:n
        if str2num(X{i}(1:2)) == curr_location
            yval = find(vec == DistributionTransformers.MAKER(i));
            xval = DistributionTransformers.MANUFDTE(i);
            plot(xval+randn(1,1)*0.1e4,yval+randn(1,1)*0.1,'k.'); hold on;
        end
    end
    ylim([0 76]);
end

figure(3);
clf;
for i = 1:length(X)
    yval = find(makers == DistributionTransformers.MAKER(i));
    if ~isempty(yval)
        xval = DistributionTransformers.MANUFDTE(i);
        if xval > 20000000
            plot((xval+randn(1,1)*0.1e4)/1e4,yval+randn(1,1)*0.1,'k.'); hold on;
        end
    end
end
xlim([1998 2021])
set(gca,'YTick',1:32,'FontSize',8)
set(gca,'XTick',2000:1:2020,'FontSize',8)
grid on
xlabel('Time','FontSize',18);
ylabel('Manufacturer','FontSize',18);
savePics('ManvsTime.pdf',1,'beamer')


figure(1); clf;
plot(DistributionTransformers.COPPERWT,DistributionTransformers.TOTALWT,'k.')
set(gca,'YScale','log','XScale','log');
hold on;
plot([100 1e6],[100 1e6],'r--');

figure(2); clf;
plot(DistributionTransformers.MANUFDTE,DistributionTransformers.COMMDTE,'k.')
hold on;
plot([min(DistributionTransformers.MANUFDTE) max(DistributionTransformers.MANUFDTE)],[min(DistributionTransformers.MANUFDTE) max(DistributionTransformers.MANUFDTE)],'r--');