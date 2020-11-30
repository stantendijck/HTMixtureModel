function [X, Y] = simulationfromQRmodel(n,QR_K,cal,QR,fiddle,TP_laplace,HS_laplace)
% input: cal: calsea, calswell, calall
% input: HT: HTsea, HTswell, HTall

X = zeros(n,1); Y = zeros(n,1);

% Simulate the indices of exceedances of the HT threshold
thr = Laplace_iCDF(fiddle.HT_quant);
U = rand(n,1);
Ithr = U < fiddle.HT_quant;
nbelow = sum(Ithr);
nexc = sum(~Ithr);

% Simulate from empirical model | X < thr
I = TP_laplace < thr;
population = [TP_laplace(I),HS_laplace(I)];
Ind = randsample(sum(I),nbelow,true);
X(Ithr) = population(Ind,1);
Y(Ithr) = population(Ind,2);

% Simulate from QR model(1) | X > thr
condX = Laplace_iCDF(U(~Ithr));

cps = cal.cps{QR_K}';

Ysim = zeros(size(condX));

tau = rand(nexc,1);
component = sum(cps < tau,2);
for i = 1:QR_K
    I = component == i;
    newtau = (tau(I) - cps(i))/(cps(i+1)-cps(i));
    prms = QR.p{QR_K}{i};
    alp = prms(1);
    gam = prms(2);
    bet = prms(3);
    Ysim(I) = gam + alp * condX(I) + condX(I).^bet .* interp1q(cal.output_p{QR_K}{i}',cal.output_q{QR_K}{i}',newtau);
end

X(~Ithr) = condX;
Y(~Ithr) = Ysim;
