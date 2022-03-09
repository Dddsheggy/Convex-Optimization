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


x_syms = sym(zeros(1));
for i=1:m
    cmd = sprintf('sym(''x%i'')',i);
    x_syms(i) = eval(cmd);
end
f = x_syms*A*inv(P)*A'*x_syms'/2 + x_syms*(A*inv(P)*q/2+b) + q'*inv(P)*A'*x_syms'/2 + q'*inv(P)*q/2;
x = ones(m,1);
temp = double(subs(f,x_syms,x'));
f_ax = [temp];
t_ax = [1];

d = (A*inv(P)*A')'*x_syms' + A*inv(P)*q/2+b + (q'*inv(P)*A')'/2;
d2 = A*inv(P)*A';
flag = 1;
k = 0; 
while(flag)
    d_temp = double(subs(d,x_syms,x'));   
    d2_temp = d2;   
    nt = -inv(d2_temp) * d_temp;
    lambda = sqrt(-d_temp' * nt);
    if(lambda^2 >= epsilon)
        t = 1;
        flag1 = 1;
        while(flag1)
           x_temp = x + t * nt;
           f_temp = double(subs(f,x_syms,x_temp'));
           if(f_temp > (temp - a * t * lambda^2))
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
           if(f_temp > (temp - a * t * lambda^2))
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
xlabel('k'),ylabel('log(f(vk)-d*)'),
title('log(f(vk)-d*)-k');