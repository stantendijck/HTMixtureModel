function allocate = getallocate(L, U)
%Get allocations based on likelihood

if nargin < 2
    U = rand(size(L,1),1);
end

allocate = cell(size(L,2),1);
I = sum([zeros(size(L, 1), 1), cumsum(L ./ sum(L, 2), 2)] < U, 2);
for ic = 1:size(L,2)
    allocate{ic} = (I == ic);
end
    
end