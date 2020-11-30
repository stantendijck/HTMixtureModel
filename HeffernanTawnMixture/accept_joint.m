function [curr_prm, accept_vec, newloglike] = accept_joint(HT, HTx_beta, curr_prm, allocate, proposed, keef)
if nargin < 6
    keef = true;
end

accept_vec = 0;

noc = size(curr_prm,2);

new_prm = proposed;
flag_prm = true; flag_keef = true;
for j = 1:noc
    %% Check if the conditions |alpha| <= 1, beta < 1, sigma >= 0 hold
    if proposed(1,j) < -1 || proposed(2,j) > 1 || proposed(4,j) < 0
        flag_prm = flag_prm & false;
    end
    if keef && proposed(1,j) > 1
        flag_prm = flag_prm & false;
    end
    %% Check the Keef Constraints
    if keef
        f = HeffernanTawn.KeefConstraints(HT.x,HT.y,new_prm(:,j),[0,1],HTx_beta(:,j));
        flag_keef = flag_keef & all(f);
    end
end

%% Update current parameters
if flag_prm && flag_keef
    %% Compute likelihood
    getRsds = @(XInd, XInd_beta, YInd, Prm)(YInd - Prm(1)*XInd - Prm(3)*XInd_beta)./(Prm(4)*XInd_beta);
    getNLOGL = @(XInd, XInd_beta, YInd, Prm)(log(Prm(4)*XInd_beta) + getRsds(XInd, XInd_beta, YInd, Prm).^2/2);
    
    NLogOld = 0;
    NLogNew = 0;
    for j = 1:noc
        xAlloc = HT.x(allocate{j});
        xAlloc_beta_old = exp(HT.logx(allocate{j}) .* curr_prm(2,j));
        xAlloc_beta_new = exp(HT.logx(allocate{j}) .* new_prm(2,j));
        yAlloc = HT.y(allocate{j});

        NLogOld = NLogOld + sum(getNLOGL(xAlloc, xAlloc_beta_old, yAlloc, curr_prm(:,j)));
        NLogNew = NLogNew + sum(getNLOGL(xAlloc, xAlloc_beta_new, yAlloc, new_prm(:,j)));
    end

    %% Accept or reject the proposed parameter
    accept = exp(NLogOld - NLogNew);
    if rand < accept
        curr_prm = proposed;
        accept_vec = 1;
        newloglike = NLogNew;
    else 
        newloglike = NLogOld;
    end
else
    newloglike = nan;
end

end