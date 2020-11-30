
clc;
I=ncinfo('NORA10_6400N_0796W.nc');

%I.Variables.Name %list of vairbale names!!

%I.Variables(10).Attributes.Name  "can find @unit scale" down here somewhere

Dt=ncread('NORA10_6400N_0796W.nc','dateser');
Dt=datetime(Dt,'convertfrom','datenum');
Hs=double(ncread('NORA10_6400N_0796W.nc','Hs'))./10;
Hs_Sea=double(ncread('NORA10_6400N_0796W.nc','Hs_Sea'))./10;
Hs_Swell=double(ncread('NORA10_6400N_0796W.nc','Hs_Swell'))./10;
Tp=double(ncread('NORA10_6400N_0796W.nc','Tp'))./10;
Tp_Sea=double(ncread('NORA10_6400N_0796W.nc','Tp_Sea'))./10;
Tp_Swell=double(ncread('NORA10_6400N_0796W.nc','Tp_Swell'))./10;

T2=double(ncread('NORA10_6400N_0796W.nc','T2'))./10;

Wind = double(ncread('NORA10_6400N_0796W.nc','WindSpeed_10m'))./10;

Drc=double(ncread('NORA10_6400N_0796W.nc','MeanWaveDirection'))./10;

Precip=double(ncread('NORA10_6400N_0796W.nc','Precipitation'))./10; 

figure(2); clf;
subplot(2,2,1)
plot(Tp,Hs,'.k')
hold on
grid on

subplot(2,2,2)
plot(Tp_Sea,Hs_Sea,'.k')
hold on
grid on

subplot(2,2,3)
plot(Tp_Swell,Hs_Swell,'.k')
hold on
grid on

subplot(2,2,4)
plot(T2,Hs,'.k')
hold on
grid on

%%
rng(12345);
T2_perm = T2 + rand(size(T2))*0.1 - 0.05;
Hs_perm = Hs + rand(size(Hs))*0.1 - 0.05;
Drc_perm = Drc + rand(size(Hs)) - 0.5;


Ind = Hs_perm./(T2_perm.^2) > 0.06;
% Ind = Hs_Sea > Hs_Swell;
figure(3);clf; plot(T2_perm(Ind),Hs_perm(Ind),'r.','MarkerSize',2); hold on;plot(T2_perm(~Ind),Hs_perm(~Ind),'k.','MarkerSize',2);
figure(4);clf; plot(T2_perm(Ind),Hs_perm(Ind)./T2_perm(Ind).^2,'r.','MarkerSize',2); hold on;plot(T2_perm(~Ind),Hs_perm(~Ind)./T2_perm(~Ind).^2,'k.','MarkerSize',2);

stormpeak = getStorms([Hs_perm,Drc,datenum(Dt),T2]);

HS_peak = Hs_perm(stormpeak(:,4));
T2_peak = T2_perm(stormpeak(:,4));
Drc_peak = Drc_perm(stormpeak(:,4));

seaHS_peak = HS_peak(Ind(stormpeak(:,4)));
swellHS_peak = HS_peak(~Ind(stormpeak(:,4)));
seaT2_peak = T2_peak(Ind(stormpeak(:,4)));
swellT2_peak = T2_peak(~Ind(stormpeak(:,4)));
seaDrc_peak = Drc_peak(Ind(stormpeak(:,4)));
swellDrc_peak = Drc_peak(~Ind(stormpeak(:,4)));

figure(5);clf; plot(seaT2_peak,seaHS_peak,'r.'); hold on; plot(swellT2_peak,swellHS_peak,'k.');

U = Laplace_iCDF(epit([T2_peak,HS_peak]));
U1 = Laplace_iCDF(epit([seaT2_peak,seaHS_peak]));
U2 = Laplace_iCDF(epit([swellT2_peak,swellHS_peak]));
figure(6);clf; subplot(1,3,1); plot(U(:,1),U(:,2),'k.');axis('square');
subplot(1,3,2); plot(U1(:,1),U1(:,2),'k.');axis('square');
subplot(1,3,3); plot(U2(:,1),U2(:,2),'k.');axis('square');


