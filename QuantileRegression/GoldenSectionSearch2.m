function x_min = GoldenSectionSearch2(fun, x0, x1, tol)

phi = (1 + sqrt(5)) / 2;

if isempty(x0) && isempty(x1)
    x_min = 0;
    return
end

x2 = x1 - (x1 - x0) / phi;
x3 = x0 + (x1 - x0) / phi;

f0 = fun(x0);
f1 = fun(x1);
f2 = fun(x2);
f3 = fun(x3);

while f0 < min([f1, f2, f3])
    x0 = x0 - (x1 - x0);
    
    x2 = x1 - (x1 - x0) / phi;
    x3 = x0 + (x1 - x0) / phi;
    
    f0 = fun(x0);
    f2 = fun(x2);
    f3 = fun(x3);
end

while f1 < min([f0, f2, f3])
    x1 = x1 + (x1 - x0);
    
    x2 = x1 - (x1 - x0) / phi;
    x3 = x0 + (x1 - x0) / phi;
    
    f1 = fun(x1);
    f2 = fun(x2);
    f3 = fun(x3);
end


while abs(x2 - x3) > tol
    if fun(x2) < fun(x3)
        x1 = x3;
    else
        x0 = x2;
    end
    
    x2 = x1 - (x1 - x0) / phi;
    x3 = x0 + (x1 - x0) / phi;
end

x_min = (x0 + x1)/2;


end