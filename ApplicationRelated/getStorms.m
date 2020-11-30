function stormpeak = getStorms(dat)

HS = dat(:,1);
Drc = dat(:,2);
SSn = dat(:,3);
Tp = dat(:,4);


%% Statistician Choices
dT = SSn(2)-SSn(1);
min_sep = 1/dT;


%% Because of ties, we add some noise to cover the rounding interval
rng(12345);
dHS = 0.1; dTp = 0.1;
HS = HS + rand(size(HS))*dHS - dHS/2;
Tp = Tp + rand(size(Tp))*dTp - dTp/2;


%% Separate storms by just looking at consecutive exceedances of threshold
peakpickingvariable = HS;

storm_threshold = quantile(peakpickingvariable,0.7);
I = peakpickingvariable > storm_threshold;

temp_storms = {};
time_counter = 1;
while time_counter < length(peakpickingvariable)
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
stormpeak = nan(length(storms),4);
for istorm = 1:length(storms)
    [~,I] = max(peakpickingvariable(storms{istorm}));
    Ind = storms{istorm}(I);
    
    D = Drc(Ind);
    if D < 45
        D = D + 360;
    end
    stormpeak(istorm,:) = [HS(Ind),D,Tp(Ind),Ind];
end

figure(1); clf;
plot(Tp,HS,'k.'); hold on;
plot(stormpeak(:,3),stormpeak(:,1),'r.');
xlabel('T_{p,ass}');ylabel('H_{S,peak}');
