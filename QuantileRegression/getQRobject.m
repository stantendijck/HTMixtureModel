function [QR, model, default] = getQRobject(input_file_name)

load(input_file_name);

model_name = model;
clear model;

% Model parameters
model.name = model_name;
model.tau_cp = tau_cp;
model.n = n;
model.large_n = 100000;
model.B = tot_iter;

% Inputs from statistician
model.constraint = model_constraint;
model.tau = tau;
model.max_n_cp = max_n_cp;

% defaults for saving/plotting intermediate/end results
default.PltOn = true;
default.saveOn = false;
default.purpose = 'poster';
default.noPrint = true;

%% Note data does not always contain bootstrapped data, sometimes also just resamples
QR = QuantileRegression();
QR.B = model.B;

for iB = 1:model.B
    QR.BS{iB} = QuantileRegression();
    QR.BS{iB}.x = data{iB}{1};
    QR.BS{iB}.y = data{iB}{2};
    
    curr_p1 = p1{iB}'; curr_p2 = p2{iB}';
    QR.BS{iB}.p = [curr_p1;curr_p2];
    
    curr_cp1 = cp1{iB}'; curr_cp2 = cp2{iB}';
    QR.BS{iB}.q = [curr_cp1;curr_cp2];
    
    QR.BS{iB}.tau = tau;
    QR.BS{iB}.model = model;
    QR.BS{iB}.default = default;
end

