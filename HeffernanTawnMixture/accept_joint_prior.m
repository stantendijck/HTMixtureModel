function [curr_prm, accept_vec, newloglike] = accept_joint_prior(HT, curr_prm, allocate, proposed, funcs)

prior = funcs{1};
from_a_to_b = funcs{2};
omega = funcs{3};

accept_vec = 0;

noc = size(curr_prm,2);
if noc > 1
	error('accept_joint_prior only supports: noc = 1');
end

%% Compute likelihood
getRsds = @(XInd, XInd_beta, YInd, Prm)(YInd - Prm(1)*XInd - Prm(3)*XInd_beta)./(Prm(4)*XInd_beta);
getNLOGL = @(XInd, XInd_beta, YInd, Prm)(log(Prm(4)*XInd_beta) + getRsds(XInd, XInd_beta, YInd, Prm).^2/2);


xAlloc = HT.x(allocate{1});
xAlloc_beta_old = exp(HT.logx(allocate{1}) .* curr_prm(2,1));
xAlloc_beta_new = exp(HT.logx(allocate{1}) .* proposed(2,1));
yAlloc = HT.y(allocate{1});

NLogOld = sum(getNLOGL(xAlloc, xAlloc_beta_old, yAlloc, curr_prm(:,1)));
NLogNew = sum(getNLOGL(xAlloc, xAlloc_beta_new, yAlloc, proposed(:,1)));


%% Accept or reject the proposed parameter
likelihoodRatio = exp(NLogOld - NLogNew);

if omega ~= 0 && omega ~= 1
	priorRatio = prior(proposed(1))/prior(curr_prm(1));
else
	if (omega == 0 && proposed(1) == 1) || (omega == 1 && proposed(1) ~= 1)
		priorRatio = 0;
	else
		priorRatio = 1;
	end
end

propNew = from_a_to_b(curr_prm(1), proposed(1));
propOld = from_a_to_b(proposed(1), curr_prm(1));
proposalRatio = propOld/propNew;

accept = likelihoodRatio * priorRatio * proposalRatio;

if rand < accept
	curr_prm = proposed;
	accept_vec = 1;
	newloglike = NLogNew;
else 
	newloglike = NLogOld;
end


end