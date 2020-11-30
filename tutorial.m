%In this tutorial, I'll show you how to use my code:
% 1. We simulate some data
% 2. We fit the QR(k) model for all 1 <= k <= 4 (approx)
% 3. We fit the HTM(k) model for input 1 <= k <= 4
% 4. How to use the fits to simulate from the models

%% First add all folders to your path
% Determine where your m-file's folder is.
folder = fileparts(which('tutorial.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% Simulate data | x > F_Laplace^{-1}(0.9)
% pick from the following:
% 1. 'log'                  simulates from evd with logistic dependenec
% 2. 'alog'                 simulates from evd with asymmetric logistic dependenec
% 3. 'normal'               simulates from gaussian copula
% 4. 'normal_log_mixture'   simulates from Distribution (F)

model.name = 'alog';
model.tau_cp = 0.5;         % mixture_probability p1 if applicable
model.n = 500;              % number of observations above thresholdp = 0.5;                    


[x, y] = DataSml2(model.name, model.n, model.tau_cp);
figure(1); clf; plot(x,y,'k.');
%% Extra defaults for fitting models

% parameters for modelling QR-model
QRmodel.tau = 0.05:0.05:0.95;           % tau_1, ..., tau_m
QRmodel.max_noc = 4;                    % max number of components to be fitted
QRmodel.B = 1;                          % number of repetitions
QRmodel.constraint = {'size',2};        % constraint that pk >= 2 * 0.05
                                        % maximum number of components - 1;
QRmodel.max_n_cp = floor(length(QRmodel.tau)/QRmodel.constraint{2}) - 1;   
              

% defaults for saving/plotting intermediate/end results
default.PltOn = true;
default.saveOn = false;
default.purpose = 'poster';
default.noPrint = true;

% HT-modelling non-exceedance probability
nep = 0.8;
u = quantile(x,nep);

%% Fit quantile-regression models

% allocation
QR = QuantileRegression(QRmodel, default, false);
QR.x = x(x>u);
QR.y = y(x>u);
QR.B = QRmodel.B;

% fit the model
fprintf('Fitting QR models\n');
[p, f, ~, q] = Fit_Model(QR.x, QR.y, QRmodel.tau, QRmodel.max_n_cp, QRmodel.constraint, default.noPrint);

QR.p = p;
QR.f = f;
QR.q = q;

% plot QR model fit
noc = 2;                % fit to be shown for QR(noc)
legendOn = false;       % show legend in fit
Plot(QR, noc, legendOn)

%% Quantile calibration stage

% specify maximum number of components that we want to use
max_noc = 4;  
pltOn = false;   % 'true': plot quantile calibration fits

% pre allocation
output_q = cell(1, max_noc);
output_p = cell(1, max_noc);
cps = cell(1, max_noc);

fprintf('Calibrating QR models\n');
for j_noc = 1:max_noc
    [output_q{j_noc}, output_p{j_noc}, cps{j_noc}] = Calibrate_Data(QR, j_noc, pltOn);
end
cal.output_q = output_q;
cal.output_p = output_p;
cal.cps = cps;

%% Fit HTM(k) model

% pre allocation
HT = cell(max_noc, 1);

% maximum number of adaptive MCMC iterations
max_iter = 10000;
for j_noc = 1:max_noc
    % Use the same data
    HT{j_noc}.x = QR.x;
    HT{j_noc}.y = QR.y;
    
    % Fit HT mixture model
    HT{j_noc} = fitHT(HT{j_noc}, max_iter, j_noc - 1, true);
end

% plot the MCMC chain for HT(noc) model
noc = 2;
plotHT(HT{noc});

%% Using the models by simulating from the models

% Simulate from QR(QR_K) model
n = 1e4;                % simulate n observations from       
QR_K = 2;               % HTM(HT_K) model     
fiddle.HT_quant = 0.8;  % with nep = fiddle.HT_quant
[X, Y] = simulationfromQRmodel(n, QR_K, cal, QR, fiddle, x, y);

figure(10); clf;
plot(X,Y,'k.');

% Simulate from HTM model
n = 1e4;                % simulate n observations from       
HT_K = 2;               % HTM(HT_K) model     
fiddle.HT_quant = 0.8;  % with nep = fiddle.HT_quant
[X, Y] = simulationfromHTmodel(n, HT_K, HT, fiddle, x, y);

figure(11); clf;
plot(X,Y,'k.');
