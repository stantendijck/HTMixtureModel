function x_min = GoldenSectionSearch(fun, a, b, tol)

if isempty(a) && isempty(b)
    x_min = 0;
    return
end


h = b - a;
invphi = (sqrt(5) - 1)/2;
invphi2 = (3 - sqrt(5))/2;

c = a + invphi2 * h;
d = a + invphi * h;

ya = fun(a);
yb = fun(b);
yc = fun(c);
yd = fun(d);

while ya < min([yb, yc, yd])
    a = a - (b - a);
    
    c = b - (b - a) * invphi;
    d = a + (b - a) * invphi;
    
    ya = fun(a);
    yc = fun(c);
    yd = fun(d);
end

while yb < min([ya, yc, yd])
    b = b + (b - a);
    
    c = b - (b - a) * invphi;
    d = a + (b - a) * invphi;
    
    yb = fun(b);
    yc = fun(c);
    yd = fun(d);
end

h = b - a;
% n = int(math.ceil(math.log(tol / h) / math.log(invphi)))

while abs(c - d) > tol
    if yc < yd
        b = d;
        d = c;
        yd = yc;
        h = invphi * h;
        c = a + invphi2 * h;
        yc = fun(c);
    else
        a = c;
        c = d;
        yc = yd;
        h = invphi * h;
        d = a + invphi * h;
        yd = fun(d);
    end
end

f = zeros(4,1);
f(1) = fun(a); 
f(2) = yc;
f(3) = yd;
f(4) = fun(b);
[~,I] = min(f);

switch I
    case 1
        x_min = a;
    case 2
        x_min = c;
    case 3
        x_min = d;
    case 4
        x_min = b;
end


end