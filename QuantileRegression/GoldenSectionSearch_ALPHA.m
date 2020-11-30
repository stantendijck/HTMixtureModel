function x_min = GoldenSectionSearch_ALPHA(fun, a, b, tol)

% phi = (1 + sqrt(5)) / 2;
% 
% if isempty(x0) && isempty(x1)
%     x_min = 0;
%     return
% end
% 
% x2 = x1 - (x1 - x0) / phi;
% x3 = x0 + (x1 - x0) / phi;
% 
% f0 = fun(x0);
% f1 = fun(x1);
% f2 = fun(x2);
% f3 = fun(x3);
% 
% while x0 > -.999 && f0 < min([f1, f2, f3])
%     x0 = max(-1, x0 - (x1 - x0));
%     
%     x2 = x1 - (x1 - x0) / phi;
%     x3 = x0 + (x1 - x0) / phi;
%     
%     f0 = fun(x0);
%     f2 = fun(x2);
%     f3 = fun(x3);
% end
% 
% while x1 < 0.999 && f1 < min([f0, f2, f3])
%     x1 = min(1, x1 + (x1 - x0));
%     
%     x2 = x1 - (x1 - x0) / phi;
%     x3 = x0 + (x1 - x0) / phi;
%     
%     f1 = fun(x1);
%     f2 = fun(x2);
%     f3 = fun(x3);
% end
% 
% 
% while abs(x2 - x3) > tol
%     if fun(x2) < fun(x3)
%         x1 = x3;
%     else
%         x0 = x2;
%     end
%     
%     x2 = x1 - (x1 - x0) / phi;
%     x3 = x0 + (x1 - x0) / phi;
% end
% 
% f0 = fun(x0);
% f1 = fun(x1);
% f2 = fun(x2);
% f3 = fun(x3);

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

while a > -.999 && ya < min([yb, yc, yd])
    a = max(-1, a - (b - a));
    
    c = b - (b - a) * invphi;
    d = a + (b - a) * invphi;
    
    ya = fun(a);
    yc = fun(c);
    yd = fun(d);
end

while b < .999 && yb < min([ya, yc, yd])
    b = min(1, b + (b - a));
    
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