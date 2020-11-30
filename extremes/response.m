function R = response(HS,TP,version)

switch version
    case 1
        alpha = 2;
        beta = 0.007;
        TP0 = 7;
    case 2
        alpha = 2;
        beta = 0.005;
        TP0 = 26;
end

R = alpha * HS ./ (1 + beta * (TP - TP0) .^ 2);

end