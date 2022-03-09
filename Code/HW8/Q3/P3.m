clc
clear all
close all
load A.mat
load b.mat
load P.mat
load q.mat

epsilon = 1e-5;
epsilon1 = 1e-14;
a = 0.2;
c = 0.8;
[m, n] = size(A);
x = A \ b;

x_syms = sym(zeros(1));
for i=1:n
    cmd = sprintf('sym(''x%i'')',i);
    x_syms(i) = eval(cmd);
end
f = x_syms * P * x_syms' / 2 + q' * x_syms';
temp = double(subs(f,x_syms,x'));
f_ax = [temp];
t_ax = [1];

d = P' * x_syms' + q;
d2 = P;
flag = 1;
k = 0; 
while(flag)
    d_temp = double(subs(d,x_syms,x'));   
    d2_temp = d2;
    nt = -inv([d2_temp A';A zeros(m,m)]) * [d_temp; zeros(m,1)];
    nt = nt(1:n);
    lambda = sqrt(nt' * d2_temp * nt);
    if(lambda^2 >= epsilon)
        t = 1;
        flag1 = 1;
        while(flag1)
           x_temp = x + t * nt;
           f_temp = double(subs(f,x_syms,x_temp'));
           if(f_temp > (temp + a * t * d_temp' * nt))
               t = c * t;
           else
               flag1 = 0;
           end
        end
        t_ax = [t_ax, t];
        x = x + t * nt;
        temp = double(subs(f,x_syms,x'));
        f_ax = [f_ax, temp];
        k = k + 1;
    elseif(lambda^2 >= epsilon1)
        t = 1;
        flag1 = 1;
        while(flag1)
           x_temp = x + t * nt;
           f_temp = double(subs(f,x_syms,x_temp'));
           if(f_temp > (temp + a * t * d_temp' * nt))
               t = c * t;
           else
               flag1 = 0;
           end
        end
        x = x + t * nt;
        temp = double(subs(f,x_syms,x'));
    else
        flag=0;
    end
end

figure(1),
plot(t_ax),
xlabel('k'),ylabel('tk'),
title('tk-k');
figure(2),
plot((1:length(f_ax)),log(f_ax - temp));
xlabel('k'),ylabel('log(f(xk)-p*)'),
title('log(f(xk)-p*)-k');