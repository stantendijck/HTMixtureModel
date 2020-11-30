function x_min = GoldenSectionSearch_BETA(fun, a, b, noPrint, tol)

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
% while abs(x2 - x3) > tol
%     if ~noPrint
%         fprintf('Curr diff = %.3f; tol = %.3f\n',abs(x2-x3),tol);
%     end
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
% f = [fun(x0),fun(x2),fun(x3),fun(x1)];
% [~,I] = min(f);
% switch I
%     case 1
%         x_min = x0;
%     case 2
%         x_min = x2;
%     case 3
%         x_min = x3;
%     case 4
%         x_min = x1;
% end

if isempty(a) && isempty(b)
    x_min = 0;
    return
end


h = b - a;
invphi = (sqrt(5) - 1)/2;
invphi2 = (3 - sqrt(5))/2;

c = a + invphi2 * h;
d = a + invphi * h;

yc = fun(c);
yd = fun(d);


% n = int(math.ceil(math.log(tol / h) / math.log(invphi)))

while abs(c - d) > tol
    if ~noPrint
        fprintf('Curr diff = %.3f; tol = %.3f\n',abs(c-d),tol);
    end

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