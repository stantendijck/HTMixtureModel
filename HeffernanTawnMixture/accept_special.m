function [curr_prm, HT] = accept_special(iter, HT, HTx_beta, curr_prm, allocate, proposed, keef)
% scheme = 1: update all separately
% scheme = 2: update all together

if nargin < 8
    keef = true;
end


[curr_prm, accept_vec] = accept_joint_special(HT, HTx_beta, curr_prm, allocate{iter}, proposed, keef);
HT.accept_vec(iter) = accept_vec;
if iter == 1
    HT.accept_rate(iter) = accept_vec;
elseif iter <= 100
    HT.accept_rate(iter) = ((iter - 1) * HT.accept_rate(iter-1) + accept_vec)/iter;
else
    HT.accept_rate(iter) = HT.accept_rate(iter-1) + (accept_vec - HT.accept_vec(iter-100))/100;
end


end