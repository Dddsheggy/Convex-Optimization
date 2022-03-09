clc
clear all
close all

syms x1 x2;
f = x1 ^ 2 / 2 + 50 * x2 ^ 2;
x0 = [100; 1];
x_hb = x0;
x_gd = x0;
epsilon = 1e-8;
a = 4/121;
b = 81/121;

x1_ax = [x_hb(1)];
x2_ax = [x_hb(2)];
f_ax = [double(subs(subs(f,x1,x_hb(1)),x2,x_hb(2)))];
f_gd_ax = f_ax;
d = -[diff(f,x1);diff(f,x2)];
flag = 1;
k = 1;

% Heavy ball 
while(flag)
    d_temp = double(subs(subs(d,x1,x_hb(1)),x2,x_hb(2)));
    nor = norm(d_temp);
    if(nor >= epsilon)
        if(k==1)
            x_hb = x_hb;
        else
            x_hb = x_hb + a * d_temp + b * (x_hb - [x1_ax(k-1); x2_ax(k-1)]);
        end
        x1_ax = [x1_ax, x_hb(1)];
        x2_ax = [x2_ax, x_hb(2)];
        f_ax = [f_ax, double(subs(subs(f,x1,x_hb(1)),x2,x_hb(2)))];
        k = k + 1;
    else
        flag=0;
    end
end
figure(1),
plot(x1_ax,x2_ax),
xlabel('x1'),ylabel('x2'),
title('x2-x1');
figure(2),
plot(log(f_ax)),
xlabel('k'),ylabel('log(f(xk))'),
title('log(f(xk))-k');

% Gradient descent 
flag = 1;
while(flag)
    t = 1;
    d_temp = double(subs(subs(d,x1,x_gd(1)),x2,x_gd(2)));
    f_org = double(subs(subs(f,x1,x_gd(1)),x2,x_gd(2)));
    nor = norm(d_temp);
    if(nor >= epsilon)
        flag1 = 1;
        while(flag1)
           x_temp = x_gd + t * d_temp;
           f_temp = double(subs(subs(f,x1,x_temp(1)),x2,x_temp(2)));
           if(f_temp > (f_org - a * t * nor ^ 2))
               t = b * t;
           else
               flag1 = 0;
           end
        end
        x_gd = x_gd + t * d_temp;
        f_gd_ax = [f_gd_ax, double(subs(subs(f,x1,x_gd(1)),x2,x_gd(2)))];
    else
        flag=0;
    end
end
figure(3),
plot(log(f_ax)),
hold on
plot(log(f_gd_ax)),
hold off
xlabel('k'),ylabel('log(f(xk))'),
title('log(f(xk))-k'),
legend('Heavy ball','Gradient descent');