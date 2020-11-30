function [curr_prm, HT] = accept(scheme, iter, HT, HTx_beta, curr_prm, allocate, proposed, keef)
% scheme = 1: update all separately
% scheme = 2: update all together

if nargin < 8
    keef = true;
end

switch scheme
    case 1
        [curr_prm, accept_vec] = accept_separate(HT, HTx_beta, curr_prm, allocate{iter}, proposed);
        HT.accept_vec{iter} = accept_vec;
    case 2
        [curr_prm, accept_vec, newlike] = accept_joint(HT, HTx_beta, curr_prm, allocate{iter}, proposed, keef);
        HT.accept_vec(iter) = accept_vec;
        if ~isnan(newlike)
            HT.Like(iter) = newlike;
        elseif iter == 1
            HT.Like(iter) = -inf;
        else
            HT.Like(iter) = HT.Like(iter-1);
        end
        if iter == 1
            HT.accept_rate(iter) = accept_vec;
        elseif iter <= 100
            HT.accept_rate(iter) = ((iter - 1) * HT.accept_rate(iter-1) + accept_vec)/iter;
        else
            HT.accept_rate(iter) = HT.accept_rate(iter-1) + (accept_vec - HT.accept_vec(iter-100))/100;
        end
    case 3
		funcs = keef;
        [curr_prm, accept_vec, newlike] = accept_joint_prior(HT, curr_prm, allocate{iter}, proposed, funcs);
        HT.accept_vec(iter) = accept_vec;
        if ~isnan(newlike)
            HT.Like(iter) = newlike;
        elseif iter == 1
            HT.Like(iter) = -inf;
        else
            HT.Like(iter) = HT.Like(iter-1);
        end
        if iter == 1
            HT.accept_rate(iter) = accept_vec;
        elseif iter <= 100
            HT.accept_rate(iter) = ((iter - 1) * HT.accept_rate(iter-1) + accept_vec)/iter;
        else
            HT.accept_rate(iter) = HT.accept_rate(iter-1) + (accept_vec - HT.accept_vec(iter-100))/100;
        end
end

end