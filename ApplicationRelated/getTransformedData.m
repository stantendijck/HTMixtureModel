function [x,y,GPD_prm,orgHS,orgTp] = getTransformedData(dat, qgpd)
% qgpd = 0.8; default
if nargin < 2
    qgpd = 0.8;
end

HS = dat(:,1);
Drc = dat(:,2);
SSn = dat(:,3);
Tp = dat(:,4);


%% Statistician Choices
dT = SSn(2)-SSn(1);
min_sep = 0.3/dT;

xl = 215; xu = 370;

q_GPD = 0.7:0.005:0.95;
[~,i_qgpd] = min(abs(q_GPD-qgpd));


%% Because of ties, we add some noise to cover the rounding interval
dHS = 0.1; dTp = 0.1;
HS = HS + rand(size(HS))*dHS - dHS/2;
Tp = Tp + rand(size(Tp))*dTp - dTp/2;


%% Separate storms by just looking at consecutive exceedances of threshold

storm_threshold = quantile(HS,0.7);
I = HS > storm_threshold;

temp_storms = {};
time_counter = 1;
while time_counter < length(HS)
    if I(time_counter)
        if ~I(time_counter-1)
            temp_storms{end+1} = time_counter;
        else
            temp_storms{end} = [temp_storms{end},time_counter];
        end
    end
    time_counter = time_counter + 1;
end

%% Merge storms if the separation is less than or equal to 8 hours

storms = {};
storms{1} = temp_storms{1};
temp_storm_counter = 2; storm_counter = 1;
while temp_storm_counter <= length(temp_storms)
    if temp_storms{temp_storm_counter}(1) - temp_storms{temp_storm_counter-1}(end) < min_sep
        storms{storm_counter} = [storms{storm_counter},(storms{storm_counter}(end)+1):(temp_storms{temp_storm_counter}(end)-1),temp_storms{temp_storm_counter}];
    else
        storm_counter = storm_counter + 1;
        storms{storm_counter} = temp_storms{temp_storm_counter};
    end
    temp_storm_counter = temp_storm_counter + 1;
end

%% Get storm peaks
stormpeak = nan(length(storms),3);
for istorm = 1:length(storms)
    [~,I] = max(HS(storms{istorm}));
    Ind = storms{istorm}(I);
    
    D = Drc(Ind);
    if D < 45
        D = D + 360;
    end
    stormpeak(istorm,:) = [HS(Ind),D,Tp(Ind)];
end

figure(1); clf;
plot(stormpeak(:,3),stormpeak(:,1),'k.');
xlabel('T_{peak}');ylabel('H_{S,peak}');

figure(2); clf;
plot(stormpeak(:,2),stormpeak(:,1),'k.');
xlabel('DRC_{peak}');ylabel('H_{S,peak}');
hold on;
plot([xl,xl],[0,14]);
plot([xu,xu],[0,14]);

%% PIT transform the observations between xl and xu in Drc
I = stormpeak(:,2) > xl & stormpeak(:,2) < xu;
data = stormpeak(I,:);

u_HS = zeros(length(data(:,1)),1);
u_Tp = zeros(length(data(:,3)),1);

GPD_threshold_HS = zeros(length(q_GPD),1);
GPD_threshold_Tp = zeros(length(q_GPD),1);
phat_HS = zeros(length(q_GPD),2);
phat_Tp = zeros(length(q_GPD),2);
for iq = 1:length(q_GPD)
    q = q_GPD(iq);
    GPD_threshold_HS(iq) = quantile(data(:,1),q);
    GPD_threshold_Tp(iq) = quantile(data(:,3),q);

    I_HS = data(:,1) > GPD_threshold_HS(iq);
    I_Tp = data(:,3) > GPD_threshold_Tp(iq);
    
    phat_HS(iq,:) = fminsearch(@(p)(GPD_like(data(I_HS,1), GPD_threshold_HS(iq),p(1),p(2))),[1,-0.1]);
    phat_Tp(iq,:) = fminsearch(@(p)(GPD_like(data(I_Tp,3), GPD_threshold_Tp(iq),p(1),p(2))),[1,-0.1]);
end

figure(3); clf;
subplot(2,2,1);
plot(q_GPD,phat_HS(:,1)); ylabel('\sigma_{HS}'); xlabel('GPD threshold');
subplot(2,2,2);
plot(q_GPD,phat_HS(:,2)); ylabel('\xi_{HS}'); xlabel('GPD threshold');
subplot(2,2,3);
plot(q_GPD,phat_Tp(:,1)); ylabel('\sigma_{Tp}'); xlabel('GPD threshold');
subplot(2,2,4);
plot(q_GPD,phat_Tp(:,2)); ylabel('\xi_{Tp}'); xlabel('GPD threshold');

q = q_GPD(i_qgpd);
I_HS = data(:,1) > GPD_threshold_HS(i_qgpd);
I_Tp = data(:,3) > GPD_threshold_Tp(i_qgpd);

prm_HS = [GPD_threshold_HS(i_qgpd),phat_HS(i_qgpd,1),phat_HS(i_qgpd,2),q];
prm_Tp = [GPD_threshold_Tp(i_qgpd),phat_Tp(i_qgpd,1),phat_Tp(i_qgpd,2),q];

u_HS(~I_HS) = epit(data(~I_HS,1))*q;
u_HS(I_HS) = GPD_CDF(data(I_HS,1),prm_HS(1),prm_HS(2),prm_HS(3))*(1-q) + q;

u_Tp(~I_Tp) = epit(data(~I_Tp,3))*q;
u_Tp(I_Tp) = GPD_CDF(data(I_Tp,3),prm_Tp(1),prm_Tp(2),prm_Tp(3))*(1-q) + q;

Laplace_HS = Laplace_iCDF(u_HS);
Laplace_Tp = Laplace_iCDF(u_Tp);

figure(4);clf;
plot(Laplace_Tp,Laplace_HS,'k.');
figure(5);clf;
subplot(1,2,1);
histogram(Laplace_Tp); xlabel('T_p');
subplot(1,2,2);
histogram(Laplace_HS); xlabel('H_S');

x = Laplace_Tp;
y = Laplace_HS;
GPD_prm = [prm_Tp;prm_HS];
orgHS = stormpeak(:,1);
orgTp = stormpeak(:,3);


end







