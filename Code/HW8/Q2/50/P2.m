clc
clear all
close all
load  A_50.mat
load d.mat
load d2.mat

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
temp = subs(f,x_syms,x');
f_ax = [double(temp)];
t_ax = [];

% 运行时间较长，存成.mat
% d = gradient(f,x_syms);
% save d
% d2 = jacobian(d,x_syms);
% save d2


flag = 1;
k = 0; 
while(flag)
    d_temp = double(subs(d,x_syms,x'));   
    d2_temp = double(subs(d2,x_syms,x'));   
    nt = -inv(d2_temp) * d_temp;
    lambda = sqrt(-d_temp' * nt);
    t = 1 / (1 + lambda);
    t_ax = [t_ax, t];
    if(lambda^2 >= epsilon)
        x = x + t * nt;
        f_temp = subs(f,x_syms,x');
        f_ax = [f_ax, double(f_temp)];
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
