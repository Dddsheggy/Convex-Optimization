clc
clear all
close all

A = load("data\1A.csv");
b = load("data\1b.csv");
[~, n] = size(A);
x = zeros(n, 1);
x_list = [x];
f_list = [f(A, b, x)];

alpha = 100;
epsilon = 2*1e-2;
epsilon_star = 1e-2;



while f_list(end) > epsilon
    x = x + inv(A'*A + eye(n)/alpha)*A'*(b - A*x);
    x_list = [x_list x];
    f_list = [f_list f(A, b, x)];
end

while f(A, b, x) > epsilon_star
    x = x + inv(A + eye(n)/alpha)*(b - A*x);
end

figure(1),
plot(log(1:length(f_list)), log(sqrt(sum((x_list - x).^2)))),
title('log(norm(x_k - x^*))-log(k)'),
xlabel('log(k)'),
ylabel('log(norm(x_k - x^*))');
figure(2),
plot(log(1:length(f_list)), log(f_list)),
title('log(f(x_k))-log(k)'),
xlabel('log(k)'),
ylabel('log(f(x_k))');
%%
function result = f(A, b, x)
    result = norm(A*x - b)^2 / 2;
end