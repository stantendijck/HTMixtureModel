function u = epit(x, y)
% Transform data x (n x d) to uniform margins
% if input y is given, then transform data to uniform margins using y as
% reference

if nargin < 2
    y = x;
end


d = size(x,2);
for i = 1:d
    x(isnan(x(:,i)),:) = [];
    y(isnan(y(:,i)),:) = [];
end

n = size(x,1); m = size(y,1);
for i = 1:d
    if length(unique(x(:,i))) ~= n
        x(:,i) = x(:,i) + (rand(n,1)-0.5)*min(diff(unique(x(:,i))));
    end
    if length(unique(y(:,i))) ~= m
        y(:,i) = y(:,i) + (rand(n,1)-0.5)*min(diff(unique(y(:,i))));
    end
end


u = repmat((1:n)',1,d);

for i = 1:d
    if isequal(x(:,i),y(:,i))
        [~,I] = sort(x(:,i));
        u(I,i) = u(:,i) / (n+1);
    else
        u(:,i) = sum(y(:,i) <= x(:,i)') / (m+1);
    end
end

end