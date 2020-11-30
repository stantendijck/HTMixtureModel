function [curr_prm, accept_vec] = accept_separate(HT, HTx_beta, curr_prm, allocate, proposed)

noc = size(curr_prm,2);
accept_vec = zeros(4,noc);

%% Loop over components
for j = 1:noc
    %% Get new proposed parameters
    for iprm = 1:4
        new_prm = curr_prm;
        new_prm(iprm,j) = proposed(iprm,j);
        %             new_prm(:,j) = proposed(:,j);
        
        %% Check if the conditions |alpha| <= 1, beta < 1, sigma >= 0 hold
        if (iprm == 1 && abs(proposed(1,j)) > 1) || (iprm == 2 && proposed(2,j) > 1) || (iprm == 4 && proposed(4,j) < 0)
            flag_prm = false;
        else
            flag_prm = true;
        end
        %% Check the Keef Constraints
        if iprm <= 2
            f = HeffernanTawn.KeefConstraints(HT.x,HT.y,new_prm(:,j),[0,1],HTx_beta(:,j));
            flag_keef = all(f);
        else
            flag_keef = true;
        end
        
        %% Update current parameters
        if flag_prm && flag_keef
            %% Compute likelihood            
            getRsds = @(XInd, XInd_beta, YInd, Prm)(YInd - Prm(1)*XInd - Prm(3)*XInd_beta)./(Prm(4)*XInd_beta);
            getNLOGL = @(XInd, XInd_beta, YInd, Prm)(log(Prm(4)*XInd_beta) + getRsds(XInd, XInd_beta, YInd, Prm).^2/2);
            
            % xAlloc = HT.x(allocate{iter, j});
            % xAlloc_beta_old = HT.x(allocate{iter, j}) .^ curr_prm(2,j);
            % xAlloc_beta_new = HT.x(allocate{iter, j}) .^ new_prm(2,j);
            % yAlloc = HT.y(allocate{iter, j});
            
            xAlloc = HT.x(allocate{j});
            xAlloc_beta_old = exp(HT.logx(allocate{j}) .* curr_prm(2,j));
            xAlloc_beta_new = exp(HT.logx(allocate{j}) .* new_prm(2,j));
            yAlloc = HT.y(allocate{j});
            
            NLogOld = getNLOGL(xAlloc, xAlloc_beta_old, yAlloc, curr_prm(:,j));
            NLogNew = getNLOGL(xAlloc, xAlloc_beta_new, yAlloc, new_prm(:,j));
            
            %% Accept or reject the proposed parameter
            accept = exp(sum(NLogOld) - sum(NLogNew));
            if rand < accept
                curr_prm(iprm,j) = proposed(iprm,j);
                accept_vec(iprm,j) = 1;
            end
        end
    end
end

end