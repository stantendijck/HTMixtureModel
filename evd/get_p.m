function p = get_p(d,alpha)
% p{i}(j) = p(i,j) in Algorithm 1.1 from Stephenson (2003)

p = cell(d,1);
p{1} = 1;
for i = 2:d
    p{i} = ones(i,1);
    for j = 1:i-1
        p{i}(1) = p{i}(1)*(1-alpha/j);
    end
    for j = 2:(i-1)
        p{i}(j) = ((i-1-alpha*j)*p{i-1}(j) + alpha*(j-1)*p{i-1}(j-1))/(i-1);
    end
    p{i}(i) = alpha^(i-1);
end


end