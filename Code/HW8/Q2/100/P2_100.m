clc
clear all
close all
load  A_100.mat
load d_100.mat

epsilon = 1e-5;
A = A_100;
[m, n] = size(A);
x = zeros(n,1);

x_syms = sym(zeros(1));
for i=1:n
    cmd = sprintf('sym(''x%i'')',i);
    x_syms(i) = eval(cmd);
end
f = -sum(log(1 - x_syms.^2));
for i=1:m
    f = f - log(1 - A(i,:)*x_syms');
end
temp = double(subs(f,x_syms,x'));
f_ax = [temp];
t_ax = [];

% 运行时间较长，已存成.mat
% d_100 = gradient(f,x_syms);
% save d_100


flag = 1;
k = 0; 
while(flag)
    d_temp = double(subs(d_100,x_syms,x'));   
    d2_temp = zeros(n);
    for i=1:n
        for j=1:n
            for h=1:m
                d2_temp(i,j) = d2_temp(i,j) + A(h,i)*A(h,j)/(1-A(h,:)*x)^2;
            end
            if(i==j)
                d2_temp(i,j) = d2_temp(i,j) + 2*(1+x(i)^2)/(1-x(i)^2)^2;
            end
        end
    end
    nt = -inv(d2_temp) * d_temp;
    lambda = sqrt(-d_temp' * nt);
    t = 1 / (1 + lambda);
    t_ax = [t_ax, t];
    if(lambda^2 >= epsilon)
        x = x + t * nt;
        f_temp = double(subs(f,x_syms,x'));
        f_ax = [f_ax, f_temp];
        k = k + 1;
    else
        flag=0;
    end
end

figure(1),
plot(log(f_ax - f_ax(k+1))),
xlabel('k'),ylabel('log(f(xk)-p*)'),
title('log(f(xk)-p*)-k');
figure(2),
plot(t_ax),
xlabel('k'),ylabel('tk'),
title('tk-k');
