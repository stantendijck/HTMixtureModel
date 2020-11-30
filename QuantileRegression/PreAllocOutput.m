if ~exist('opt_prm','var')
    alpha = [];
    beta = [];
    z = [];
    opt_prm = [];
    fval = 0;
end

if ~exist('opt_prm_HT','var')
    alpha_HT = [];
    beta_HT = [];
    z_HT = [];
    opt_prm_HT = [];
    fval_HT = 0;
end

if ~exist('opt_prm_HT_median','var')
    alpha_HT_median = [];
    beta_HT_median = [];
    z_HT_median = [];
    opt_prm_HT_median = [];
    fval_HT_median = 0;
end

if ~exist('opt_prm_M','var')
    alpha_M = [];
    beta_M = [];
    z_M = [];
    opt_prm_M = [];
    fval_M = 0;
    opt_cp = -1;
    tau_M = [];
end

if ~exist('opt_prm_HT_M','var')
    alpha_HT_M = [];
    beta_HT_M = [];
    z_HT_M = [];
    opt_prm_HT_M = [];
    fval_HT_M = 0;
    opt_cp_HT = -1;
    tau_HT_M = [];
end
    
if ~exist('opt_prm_HT_median_M','var')
    alpha_HT_median_M = [];
    beta_HT_median_M = [];
    z_HT_median_M = [];
    opt_prm_HT_median_M = [];
    fval_HT_median_M = 0;
    opt_cp_HT_median = -1;
    tau_HT_median_M = [];
end
