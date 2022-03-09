clc
clear all
close all
load ../data/A.mat
load ../data/b.mat
load ../data/P.mat
load ../data/q.mat
load ../data/x_0.mat
load ../data/mu.mat

epsilon = 1e-8;
alpha = 0.1;
beta = 0.8;
[m, n] = size(A);
x = x_0;
v = mu;
mu = 10;
t = 10;
t_ax = [t];

while(n/t >= epsilon)
    d_temp = t*P*x + t*q - 1./abs(x);
    d2_temp = t*P + diag(1./x.^2);
    xv = -inv([d2_temp, A';A, zeros(m,m)]) * [d_temp + A'*v; A*x - b];
    nt = xv(1:n);
    dv = xv(n+1:end);
    r = [d_temp+A'*v; A*x - b];
    
    while(norm(r) >= epsilon*t)
        t_nt = 1;
        flag1 = 1;
        while(flag1)
            x_temp = x + t_nt * nt;
            v_temp = v + t_nt * dv;
            r_temp = [t*P*x_temp+t*q-1./abs(x_temp) + A'*v_temp; A*x_temp - b];
            if(norm(r_temp) > (1-alpha*t_nt)*norm(r) || sum(x_temp<0) ~= 0)
                t_nt = beta * t_nt;
            else
                flag1 = 0;
            end
        end
        x = x_temp;
        v = v_temp;
        
        d_temp = t*P*x + t*q - 1./abs(x);
        d2_temp = t*P + diag(1./x.^2);
        xv = -inv([d2_temp, A';A, zeros(m,m)]) * [d_temp + A'*v; A*x - b];
        nt = xv(1:n);
        dv = xv(n+1:end);
        r = [d_temp+A'*v; A*x - b];
        t_ax = [t_ax, t];
        fx = x'*P*x/2 + q'*x; 
    end    
    t = mu * t; 
end

% 存储结果
lambda=1./t*x;
save('res_P2_1.mat','x','v','lambda');

figure(1),
plot(log(n./t_ax)),
xlabel('k'),ylabel('log(n/t)'),
title('log(n/t)-k');