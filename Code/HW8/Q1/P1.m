clc
clear
close 

syms x1 x2;
f = (10 * x1^2 + x2^2) / 2 + 5 * log(1 + exp(-x1 - x2));
x = [0; 0];
epsilon = 1e-8;
a = 0.2;
b = 0.8;
x1_ax = [x(1)];
x2_ax = [x(2)];
f_ax = [double(subs(subs(f,x1,x(1)),x2,x(2)))];

d = [diff(f,x1);diff(f,x2)];
d2 = [diff(d,x1) diff(d,x2)]; 
flag = 1;
k = 0; 
while(flag)
    t = 1;
    d_temp = double(subs(subs(d,x1,x(1)),x2,x(2)));
    d2_temp = double(subs(subs(d2,x1,x(1)),x2,x(2)));
    nt = -inv(d2_temp) * d_temp;
    f_org = double(subs(subs(f,x1,x(1)),x2,x(2)));
    nor = double(norm(d_temp));
    if(nor >= epsilon)
        flag1 = 1;
        while(flag1)
           x_temp = x + t * nt;
           f_temp = double(subs(subs(f,x1,x_temp(1)),x2,x_temp(2)));
           if(f_temp > (f_org + a * t * d_temp' * nt))
               t = b * t;
           else
               flag1 = 0;
           end
        end
        x = x + t * nt;
        x1_ax = [x1_ax, x(1)];
        x2_ax = [x2_ax, x(2)];
        f_ax = [f_ax, double(subs(subs(f,x1,x(1)),x2,x(2)))];
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