% Determine where your folder is.
folder = fileparts(which('runApplication.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


%% Read in data
readdata;

Stormpeak = [HS_peak,Drc_peak,T2_peak];
seaStormpeak = [seaHS_peak,seaDrc_peak,seaT2_peak];
swellStormpeak = [swellHS_peak,swellDrc_peak,swellT2_peak];

%% Model - fixed parameters

% modelling parameters
model.tau = 0.05:0.05:0.95;
model.max_noc = 2;
model.B = 1;
model.constraint = {'size',2};
model.max_n_cp = 9;
model.name = 'data';

% defaults for saving/plotting intermediate/end results
default.PltOn = false;
default.saveOn = false;
default.purpose = 'poster';
default.noPrint = true;

% fiddle parameters
fiddle.B = 1; fiddle.marg_quant = 0.7; fiddle.HT_quant = 0.8;

%% The data analysis

%% Transform data to Laplace margins

% sea waves
[TP_laplace_sea, HS_laplace_sea, GPD_prm_sea, HS_original_sea, TP_original_sea] = PIT(seaStormpeak, fiddle.marg_quant);
% swell waves
[TP_laplace_swell, HS_laplace_swell, GPD_prm_swell, HS_original_swell, TP_original_swell] = PIT(swellStormpeak, fiddle.marg_quant);
% all waves simultaneously
[TP_laplace, HS_laplace, GPD_prm, HS_original, TP_original] = PIT(Stormpeak, fiddle.marg_quant);

% one component analysis on sea and swell waves
model.max_noc = 1;
fprintf('Fitting one component model for sea waves\n');
[pB_median_after_burnInSEA, estimates_matSEA, A_originalSEA, QRsea, HTsea, calsea] = DataAnalysisWrapper(seaStormpeak, fiddle, model, default);
fprintf('Fitting one component model for swell waves\n');
[pB_median_after_burnInSWELL, estimates_matSWELL, A_originalSWELL, QRswell, HTsswell, calswell] = DataAnalysisWrapper(swellStormpeak, fiddle, model, default);

% two component analysis on all waves simultaneously
model.max_noc = 2;
fprintf('Fitting two component model for all waves simultaneously\n');
[pB, est, A_orig, QRall, HTall, calall] = DataAnalysisWrapper(Stormpeak, fiddle, model, default);


%% Plot data analysis results
figure(1); clf;
plot(T2_perm,Hs_perm,'k.'); hold on;
plot(TP_original,HS_original,'r.','MarkerSize',10);
xlabel('T_{2}'); ylabel('H_{S}');

clear x y estimates_mat
x.sea = TP_original_sea; x.swell = TP_original_swell;
y.sea = HS_original_sea; y.swell = HS_original_swell;
estimates_mat.sea = estimates_matSEA; estimates_mat.swell = estimates_matSWELL;
pB_mat.sea = pB_median_after_burnInSEA; pB_mat.swell = pB_median_after_burnInSWELL;

colorbarflag = false;
figure(3); clf;
subplot(3,2,1);
plotResults(TP_original, HS_original, 1, est, A_orig, false);
xlim([6.3 16.6]);  ylim([3.6 18.1]);
title('QR(1)');
set(gca,'XTick',[],'YTick',[]);
set(gca,'FontSize',15);

subplot(3,2,2);
plotResults(TP_original, HS_original, 1, pB, A_orig, false);
xlim([6.3 16.6]);  ylim([3.6 18.1]);
title('HT(1)');
set(gca,'XTick',[],'YTick',[]);
set(gca,'FontSize',15);

subplot(3,2,3);
plotResults_PartitioningMethod(x, y, 1, estimates_mat, A_orig, false);
xlim([6.3 16.6]);  ylim([3.6 18.1]);
title('Part-QR(1)');
set(gca,'XTick',[],'YTick',[]);
set(gca,'FontSize',15);

subplot(3,2,4);
plotResults_PartitioningMethod(x, y, 1, pB_mat, A_orig, false);
xlim([6.3 16.6]);  ylim([3.6 18.1]);
title('Part-HT(1)');
set(gca,'XTick',[],'YTick',[]);
set(gca,'FontSize',15);

subplot(3,2,5);
plotResults_QR(TP_original, HS_original, 2, est, A_orig, false);
xlabel('T_{2,ass}'); ylabel('H_{S,peak}');
xlim([6.3 16.6]);  ylim([3.6 18.1]);
title('QR(2)');
set(gca,'FontSize',15);

subplot(3,2,6);
plotResults(TP_original, HS_original, 2, pB, A_orig, false);
xlim([6.3 16.6]);  ylim([3.6 18.1]);
title('HT(2)');
set(gca,'XTick',[],'YTick',[]);
set(gca,'FontSize',15);


%% Response calculations - results/plots

%% Initialisation

% Generic synthetic response function
fun = @(hs,tp,prm)(prm.alpha * hs ./ (1 + prm.beta * abs(tp - prm.tp0).^(2)));
% Response R1 in the paper
prm.alpha = 2; prm.beta = 1; prm.tp0 = 16;
% Response R2 in the paper
prm_2.alpha = 2; prm_2.beta = 2; prm_2.tp0 = 26;

%% Simulate from models (for estimation of response distributions)

n = 1e6; % number of simulated observations

% QR model
% simulate from one component models: sea waves, swell waves and all waves simultaneously
[XsimseaQR, YsimseaQR] = simulationfromQRmodel(n,1,calsea,QRsea,fiddle,TP_laplace_sea,HS_laplace_sea);
[XsimswellQR, YsimswellQR] = simulationfromQRmodel(n,1,calswell,QRswell,fiddle,TP_laplace_swell,HS_laplace_swell);
[XsimallQR, YsimallQR] = simulationfromQRmodel(n,1,calall,QRall,fiddle,TP_laplace,HS_laplace);
% simulate from two component models: all waves simultaneously
[XsimallQR2, YsimallQR2] = simulationfromQRmodel(n,2,calall,QRall,fiddle,TP_laplace,HS_laplace);

% HT model
% simulate from one component models: sea waves, swell waves and all waves simultaneously
[XsimseaHT, YsimseaHT] = simulationfromHTmodel(n,1,HTsea,fiddle,TP_laplace_sea,HS_laplace_sea);
[XsimswellHT, YsimswellHT] = simulationfromHTmodel(n,1,HTsswell,fiddle,TP_laplace_swell,HS_laplace_swell);
[XsimallHT, YsimallHT] = simulationfromHTmodel(n,1,HTall,fiddle,TP_laplace,HS_laplace);
% simulate from two component models: all waves simultaneously
[XsimallHT2, YsimallHT2] = simulationfromHTmodel(n,2,HTall,fiddle,TP_laplace,HS_laplace);


% We transform each sample to original margins
TP_original_sea_sim_QR = transform2original(XsimseaQR,GPD_prm_sea(1,:),TP_original_sea);
HS_original_sea_sim_QR = transform2original(YsimseaQR,GPD_prm_sea(2,:),HS_original_sea);

TP_original_swell_sim_QR = transform2original(XsimswellQR,GPD_prm_swell(1,:),TP_original_swell);
HS_original_swell_sim_QR = transform2original(YsimswellQR,GPD_prm_swell(2,:),HS_original_swell);

TP_original_all_sim_QR = transform2original(XsimallQR,GPD_prm(1,:),TP_original);
HS_original_all_sim_QR = transform2original(YsimallQR,GPD_prm(2,:),HS_original);

TP_original_all_sim_QR2 = transform2original(XsimallQR2,GPD_prm(1,:),TP_original);
HS_original_all_sim_QR2 = transform2original(YsimallQR2,GPD_prm(2,:),HS_original);

TP_original_sea_sim_HT = transform2original(XsimseaHT,GPD_prm_sea(1,:),TP_original_sea);
HS_original_sea_sim_HT = transform2original(YsimseaHT,GPD_prm_sea(2,:),HS_original_sea);

TP_original_swell_sim_HT = transform2original(XsimswellHT,GPD_prm_swell(1,:),TP_original_swell);
HS_original_swell_sim_HT = transform2original(YsimswellHT,GPD_prm_swell(2,:),HS_original_swell);

TP_original_all_sim_HT = transform2original(XsimallHT,GPD_prm(1,:),TP_original);
HS_original_all_sim_HT = transform2original(YsimallHT,GPD_prm(2,:),HS_original);

TP_original_all_sim_HT2 = transform2original(XsimallHT2,GPD_prm(1,:),TP_original);
HS_original_all_sim_HT2 = transform2original(YsimallHT2,GPD_prm(2,:),HS_original);


% Combine sea and swell waves together
seaprob = length(TP_original_sea)/length(TP_original);
I = rand(n,1) < seaprob;
comb_model_TP_QR = zeros(n,1); 
comb_model_HS_QR = zeros(n,1); 
comb_model_TP_HT = zeros(n,1); 
comb_model_HS_HT = zeros(n,1); 

comb_model_TP_QR(I) = TP_original_sea_sim_QR(I);
comb_model_HS_QR(I) = HS_original_sea_sim_QR(I);
comb_model_TP_QR(~I) = TP_original_swell_sim_QR(~I);
comb_model_HS_QR(~I) = HS_original_swell_sim_QR(~I);

comb_model_TP_HT(I) = TP_original_sea_sim_HT(I);
comb_model_HS_HT(I) = HS_original_sea_sim_HT(I);
comb_model_TP_HT(~I) = TP_original_swell_sim_HT(~I);
comb_model_HS_HT(~I) = HS_original_swell_sim_HT(~I);

%% synthetic response R1
data_allQR = fun(TP_original_all_sim_QR, HS_original_all_sim_QR, prm);

manual_linspace_vec = linspace(max(0,min(data_allQR)-1),max(data_allQR)*1.5,1000);

p_allQR = []; cntr = 1;
for i = manual_linspace_vec
    p_allQR(cntr) = log(mean(data_allQR>i))/log(10);
    cntr = cntr + 1;
end

data_allQR2 = fun(TP_original_all_sim_QR2, HS_original_all_sim_QR2, prm);
p_allQR2 = []; cntr = 1;
for i = manual_linspace_vec
    p_allQR2(cntr) = log(mean(data_allQR2>i))/log(10);
    cntr = cntr + 1;
end

data_allHT = fun(TP_original_all_sim_HT, HS_original_all_sim_HT, prm);
p_allHT = []; cntr = 1;
for i = manual_linspace_vec
    p_allHT(cntr) = log(mean(data_allHT>i))/log(10);
    cntr = cntr + 1;
end

data_allHT2 = fun(TP_original_all_sim_HT2, HS_original_all_sim_HT2, prm);
p_allHT2 = []; cntr = 1;
for i = manual_linspace_vec
    p_allHT2(cntr) = log(mean(data_allHT2>i))/log(10);
    cntr = cntr + 1;
end

data_QRcomb = fun(comb_model_TP_QR, comb_model_HS_QR, prm);
p_QRcomb = []; cntr = 1;
for i = manual_linspace_vec
    p_QRcomb(cntr) = log(mean(data_QRcomb>i))/log(10);
    cntr = cntr + 1;
end

data_HTcomb = fun(comb_model_TP_HT, comb_model_HS_HT, prm);
p_HTcomb = []; cntr = 1;
for i = manual_linspace_vec
    p_HTcomb(cntr) = log(mean(data_HTcomb>i))/log(10);
    cntr = cntr + 1;
end

%% synthetic response R2
data_allQR_2 = fun(TP_original_all_sim_QR, HS_original_all_sim_QR, prm_2);

manual_linspace_vec_2 = linspace(max(0,min(data_allQR_2)-1),max(data_allQR_2)*1.5,1000);

p_allQR_2 = []; cntr = 1;
for i = manual_linspace_vec_2
    p_allQR_2(cntr) = log(mean(data_allQR_2>i))/log(10);
    cntr = cntr + 1;
end

data_allQR2_2 = fun(TP_original_all_sim_QR2, HS_original_all_sim_QR2, prm_2);
p_allQR2_2 = []; cntr = 1;
for i = manual_linspace_vec_2
    p_allQR2_2(cntr) = log(mean(data_allQR2_2>i))/log(10);
    cntr = cntr + 1;
end

data_allHT_2 = fun(TP_original_all_sim_HT, HS_original_all_sim_HT, prm_2);
p_allHT_2 = []; cntr = 1;
for i = manual_linspace_vec_2
    p_allHT_2(cntr) = log(mean(data_allHT_2>i))/log(10);
    cntr = cntr + 1;
end

data_allHT2_2 = fun(TP_original_all_sim_HT2, HS_original_all_sim_HT2, prm_2);
p_allHT2_2 = []; cntr = 1;
for i = manual_linspace_vec_2
    p_allHT2_2(cntr) = log(mean(data_allHT2_2>i))/log(10);
    cntr = cntr + 1;
end

data_QRcomb_2 = fun(comb_model_TP_QR, comb_model_HS_QR, prm_2);
p_QRcomb_2 = []; cntr = 1;
for i = manual_linspace_vec_2
    p_QRcomb_2(cntr) = log(mean(data_QRcomb_2>i))/log(10);
    cntr = cntr + 1;
end

data_HTcomb_2 = fun(comb_model_TP_HT, comb_model_HS_HT, prm_2);
p_HTcomb_2 = []; cntr = 1;
for i = manual_linspace_vec_2
    p_HTcomb_2(cntr) = log(mean(data_HTcomb_2>i))/log(10);
    cntr = cntr + 1;
end

%% Transform the above to return levels

m = length(TP_original); n = 530292/(24*365.25);

% R1
ReturnPeriod_allQR =  m/n * 10.^(-p_allQR);
ReturnPeriod_allQR2 =  m/n * 10.^(-p_allQR2);
ReturnPeriod_allHT =  m/n * 10.^(-p_allHT);
ReturnPeriod_allHT2 =  m/n * 10.^(-p_allHT2);
ReturnPeriod_QRcomb =  m/n * 10.^(-p_QRcomb);
ReturnPeriod_HTcomb =  m/n * 10.^(-p_HTcomb);

% R2
ReturnPeriod_allQR_2 =  m/n * 10.^(-p_allQR_2);
ReturnPeriod_allQR2_2 =  m/n * 10.^(-p_allQR2_2);
ReturnPeriod_allHT_2 =  m/n * 10.^(-p_allHT_2);
ReturnPeriod_allHT2_2 =  m/n * 10.^(-p_allHT2_2);
ReturnPeriod_QRcomb_2 =  m/n * 10.^(-p_QRcomb_2);
ReturnPeriod_HTcomb_2 =  m/n * 10.^(-p_HTcomb_2);

%% Plot all results for R1

figure(100); clf
for k = 1:6
    linestyle = '-';
    switch k
        case 1
            subplot(3,2,k);
            plot(ReturnPeriod_allQR,manual_linspace_vec,'k','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_1'); title('QR(1)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 35]);
        case 2
            subplot(3,2,k);
            plot(ReturnPeriod_allQR2,manual_linspace_vec,'r','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_1'); title('QR(2)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 35]);
        case 3
            subplot(3,2,k);
            plot(ReturnPeriod_allHT,manual_linspace_vec,'g','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_1'); title('HT(1)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 35]);
        case 4
            subplot(3,2,k);
            plot(ReturnPeriod_allHT2,manual_linspace_vec,'b','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_1'); title('HT(2)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 35]);
        case 5
            subplot(3,2,k);
            plot(ReturnPeriod_QRcomb,manual_linspace_vec,'m','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_1'); title('QRcombined');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 35]);
        case 6
            subplot(3,2,k);
            plot(ReturnPeriod_HTcomb,manual_linspace_vec,'c','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_1'); title('HTcombined');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 35]);
    end
end

figure(101); clf
for k = 1:6
    linestyle = '-';
    switch k
        case 1
            subplot(3,2,k);
            plot(ReturnPeriod_allQR_2,manual_linspace_vec_2,'k','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_2'); title('QR(1)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 0.3]);
        case 2
            subplot(3,2,k);
            plot(ReturnPeriod_allQR2_2,manual_linspace_vec_2,'r','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_2'); title('QR(2)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 0.3]);
        case 3
            subplot(3,2,k);
            plot(ReturnPeriod_allHT_2,manual_linspace_vec_2,'g','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_2'); title('HT(1)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 0.3]);
        case 4
            subplot(3,2,k);
            plot(ReturnPeriod_allHT2_2,manual_linspace_vec_2,'b','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R'); title('HT(2)');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 0.3]);
        case 5
            subplot(3,2,k);
            plot(ReturnPeriod_QRcomb_2,manual_linspace_vec_2,'m','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_2'); title('QRcombined');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 0.3]);
        case 6
            subplot(3,2,k);
            plot(ReturnPeriod_HTcomb_2,manual_linspace_vec_2,'c','LineStyle',linestyle); hold on;
            xlabel('Return period'); ylabel('R_2'); title('HTcombined');
            set(gca,'XScale','log'); xx = xlim;
            xlim([10 1e7]); ylim([0 0.3]);
    end
end



