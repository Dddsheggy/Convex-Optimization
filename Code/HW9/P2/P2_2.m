clc
clear all
close all
load ./data/A.mat
load ./data/b.mat
load ./data/P.mat
load ./data/q.mat
load ./data/x_0.mat
load ./data/lambda.mat
load ./data/mu.mat

epsilon = 1e-8;
alpha = 0.1;
beta = 0.8;
[m, n] = size(A);
x = x_0;
Df = -eye(n);
eta_ax = [];
r_feas_ax = [];
const = 10;

flag = 1;
while(flag)
    % calculate t
    f = x'*P*x/2 + q'*x;
    eta = -(-x)' * lambda;
    eta_ax = [eta_ax, eta];
    t = const * n / eta;
    
    % calculate delta y_pd
    r_dual = P*x + q + Df'*lambda + A'*mu;
    r_cent = -diag(lambda)*(-x) - 1/t*ones(n,1);
    r_pri = A*x - b;
    r_feas = sqrt(norm(r_pri)^2 + norm(r_dual)^2);
    r_feas_ax = [r_feas_ax, r_feas];
    pd = -inv([P, Df', A'; -diag(lambda)*Df, -diag(-x), zeros(n,m); A, zeros(m,n+m)])...
        * [r_dual; r_cent; r_pri];
    x_pd = pd(1:n);
    lambda_pd = pd(n+1:2*n);
    mu_pd = pd(2*n+1:2*n+m);
    
    % update
    if(norm(r_pri)<=epsilon && norm(r_dual)<=epsilon && eta<=epsilon)
        flag = 0;
    else
        % calculate s
        temp = [];
        for i=1:n
            if(lambda_pd(i)<0)
                temp = [temp, -lambda(i)/lambda_pd(i)];
            end
        end
        s = 0.99 * min(1,min(temp));
        flag1 = 1;
        while(flag1)
            f_plus = -(x+s*x_pd);
            if(f_plus < 0)
                flag1 = 0;
            else
                s = beta * s;                             
            end
        end

        flag1 = 1;
        rt = [r_dual; r_cent; r_pri];
        while(flag1)
            rt_plus = [P*(x+s*x_pd) + q + Df'*(lambda+s*lambda_pd) + A'*(mu+s*mu_pd);...
                -diag(lambda+s*lambda_pd)*(-(x+s*x_pd)) - 1/t*ones(n,1); A*(x+s*x_pd) - b];
            if(norm(rt_plus)<=(1-alpha*s)*norm(rt))
                flag1 = 0;
            else
                s = beta * s;
            end
        end
        
        x = x + s*x_pd;
        lambda = lambda + s*lambda_pd;
        mu = mu + s*mu_pd;
    end  
end

v = mu;
% 存储结果
save('res_P2_2.mat','x','v','lambda');

figure(1),
plot(log(eta_ax)),
xlabel('k'),ylabel('log(\eta)'),
title('log(\eta)-k');
figure(2),
plot(log(r_feas_ax)),
xlabel('k'),ylabel('log(r_{feas})'),
title('log(r_{feas})-k');
