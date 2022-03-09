clc
clear all
close all

A = load("data\2A.csv");
b = load("data\2b.csv");
x = load("data\2x0.csv");
x = x';
c = 0.01;
beta = 0.5;
epsilon = 1e-8;
epsilon_star = 1e-9;
iter_norm_list = [];
x_list = [x x];


flag = 1;
k = 1; 
while flag
    former_x = x;
    alphak = c * power(k, -beta);
    test = gk(A, b, x) - A'*(A*x-b);
    x = x - alphak*gk(A, b, x);
    iter_norm = norm(x - former_x);
    iter_norm_list = [iter_norm_list iter_norm];
    x_list = [x_list x];
    k = k + 1;
    if(iter_norm^2 < epsilon)
        flag = 0;
    end
end

h_list = [];
for i=1:size(x_list,2)
    h_list = [h_list h(A, b, x_list(:,i))];
end


while iter_norm^2 > epsilon_star
    former_x = x;
    alphak = c * power(k, -beta);
    x = x - alphak*gk(A, b, x);
    iter_norm = norm(x - former_x);
    k = k + 1;
end

figure(1),
plot(log(1:length(iter_norm_list)), log(iter_norm_list)),
title('log(norm(x_{k+1} - x_k))-log(k)'),
xlabel('log(k)'),
ylabel('log(norm(x_{k+1} - x_k))');
figure(2),
plot(log(1:length(h_list)), log(sqrt(sum((x_list - x).^2)))),
title('log(norm(x_k - x^*))-log(k)'),
xlabel('log(k)'),
ylabel('log(norm(x_k - x^*))');
figure(3),
plot(log(1:length(h_list)), log(h_list)),
title('log(h(x_k))-log(k)'),
xlabel('log(k)'),
ylabel('log(h(x_k))');
%%
function result = h(A, b, x)
    result = norm(A*x - b)^2/2 + norm(x, 1);
end
function result = gk(A, b, x)
    temp = x;
    temp(x>0) = 1;
    temp(x<0) = -1;
    temp(abs(x)<1e-5) = 1 - 2*rand(length(find(abs(x)<1e-5)), 1);
    result = A'*(A*x-b) + temp;
end