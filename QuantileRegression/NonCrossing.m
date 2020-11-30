function flag = NonCrossing(param, u, tau)
%This function checks whether the parametric quantiles stored in param are
%crossing. If they do not cross, this function returns flag = True
%
% Example:
% a = [0.4,0.5];
% b = [1.2,0.5];
% c = [0.5,0.4];
% u = 0;
%
% x_vec = u:0.1:(10*u+10);
% 
% y_1 = (a(1)-a(2))*x_vec + b(1)*x_vec.^c(1) - b(2)*x_vec.^c(2);
% y_2 = (a(1)-a(2))*x_vec.^(1-c(2)) + b(1)*x_vec.^(c(1)-c(2)) - b(2);
% 
% figure(1); clf;
% subplot(1,2,1);
% plot(x_vec,y_1)
% subplot(1,2,2);
% plot(x_vec,y_2)

flag = 1;

if iscell(param)
    n_cp = length(param{2});
    
    t = cell(n_cp+1,1);
    if n_cp == 0
        t{1} = tau;
    else
        t{1} = tau(tau<=param{2}(1));
        for i = 1:n_cp-1
            t{i+1} =  tau(tau>param{2}(i)&tau<=param{2}(i+1));
        end
        t{end} = tau(tau>param{2}(end));
    end

    length_tau = length(tau);
    a = zeros(length_tau,1);
    b = param{1}(1+2*length(t):end);
    c = zeros(length_tau,1);
    
    cnt = 1;
    for i = 1:length(t)
        a(cnt:cnt+length(t{i})-1) = param{1}(i);
        c(cnt:cnt+length(t{i})-1) = param{1}(i+length(t));
        cnt = cnt+length(t{i});
    end
else
    a = param(1,:);
    b = param(2,:);
    c = param(3,:);
end

dim = length(a);
if dim == 1    
    return
end

i = 1;
while i <= dim - 1 && flag
    if a(i) > a(i+1)
        flag = 0;
    else
        if sign(c(i) - c(i+1)) == sign(b(i))
            x_pow = (1-c(i+1))*(a(i+1)-a(i))/(c(i)-c(i+1))/b(i);
            x_opt = x_pow^(1/(c(1)-1));
            if x_opt < u
                check = (a(i)-a(i+1))*u + b(i)*u.^c(i) - b(i+1)*u.^c(i+1);
            else
                check = (a(i+1)-a(i))*x_pow^((1-c(i+1))/(c(i)-1)) + b(i)*x_pow^((c(i)-c(i+1))/(c(i)-1)) - b(i+1);
            end
        elseif sign(c(i)-c(i+1))==0 && b(i) < b(i+1)
            check = -1;
        else
            check = (a(i)-a(i+1))*u + b(i)*u.^c(i) - b(i+1)*u.^c(i+1);
        end
        flag = check <= 0;
    end
    i = i + 1;
end

end