function proposed = propose(curr_prm, scheme, stepsize)

nop = size(curr_prm,1);
noc = size(curr_prm,2);
switch scheme
    case 1
        proposed = curr_prm + randn(nop,noc) .* stepsize.h;
    case 2
        proposed = reshape(reshape(curr_prm, nop*noc, 1)+(1-0.05)*mvnrnd(zeros(nop*noc,1),stepsize.adp_h,1)'*2.38/sqrt(nop*noc)+0.05*reshape(stepsize.h1,nop*noc,1)/sqrt(nop*noc).*randn(nop*noc,1),nop,noc);
end

end