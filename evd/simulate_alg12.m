function [x, y] = simulate_alg12(n, d, alpha, theta)
% simulate from asymmetric logistic model using Algorithm 1.2 in Stephenson
% (2003)
%
% INPUT:
% n         number of samples
% d         dimension
% alpha     (2^d - d - 1) x 1 dependence parameters
% theta     d*2^(d-1) x 1 dependence parameters
%
% OUTPUT:
% x         (n x d) sample on Gumbel scale
% y         (n x d) sample on Frechet scale

% Check input
if length(alpha) ~= 2^d - d - 1
    error('In simulate_alg12, length(alpha) needs to be equal to 2^d - d - 1');
end
if length(theta) ~= d*2^(d-1)
    error('In simulate_alg12, length(theta) needs to be equal to d*2^(d-1)');
end
if any(alpha) <=0  || any(alpha > 1)
    error('In simulate_alg12, 0 < alpha(i) <= 1 for all i');
end

B = cell(1,d);
for iB = 1:d
    B{iB} = combnk(1:d,iB);
end

for i = 1:d
    s = 0;
    counter = 1;
    for j = 1:d
        for k = 1:size(B{j},1)
            for m = 1:size(B{j},2)
                if B{j}(k,m) == i
                    s = s + theta(counter);
                end
                counter = counter + 1;
            end
        end
    end
    if s ~= 1
        error('In simulate_alg12, sum of theta_%d over all sets b does not equal 1',i);
    end
end

flag = false;
counter_alpha = 1;
counter_theta = d + 1;
for i = 2:d
    for j = 1:size(B{i},1)
        if alpha(counter_alpha) == 1
            for k = 1:size(B{i},2)
                if theta(counter_theta) ~= 0
                    flag = true;
                end
                counter_theta = counter_theta + 1;
            end
        end
        counter_alpha = counter_alpha + 1;
    end
end
if flag
    error('In simulate_alg12, if alpha(b) = 1, then theta(i,b) = 0 for all i must hold');
end

flag = false;
counter_theta = d + 1;
for i = 2:d
    for j = 1:size(B{i},1)
        s = 0;
        for k = 1:size(B{i},2)
            s = s + (theta(counter_theta)==0);
            counter_theta = counter_theta + 1;
        end
        if s == 1
            flag = true;
        end
    end
end
if flag
    error('In simulate_alg12, if theta(j,b)=0 for all j in b_{-i}, then theta(i,b) needs to be 0');
end

x_b = cell(n,d);

for i_n = 1:n
    counter_alpha = 1;
    counter_theta = 1;
    t_alpha = 1;
    for i = 1:d
        x_b{i_n,i} = zeros(size(B{i}));
        for j = 1:size(B{i},1)
            if i >= 2
                t_alpha = alpha(counter_alpha);
                counter_alpha = counter_alpha + 1;
            end
            [~,z] = simulate_alg11(1, size(B{i},2), t_alpha);
            for k = 1:size(B{i},2)
                x_b{i_n,i}(j,k) = theta(counter_theta)*z(k);
                counter_theta = counter_theta + 1;
            end
        end
    end
end

y = zeros(n,d);
for i_n = 1:n
    y(i_n,:) = -Inf;
    for i = 1:d
        for j = 1:d
            for k = 1:size(B{j},1)
                for m = 1:size(B{j},2)
                    if B{j}(k,m) == i
                        y(i_n,i) = max(y(i_n,i),x_b{i_n,j}(k,m));
                    end
                end
            end
        end
    end
end

x = log(y);


end