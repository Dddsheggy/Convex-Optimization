clc
clear all
close all

A1 = load("problem1data\A1.csv");
A2 = load("problem1data\A2.csv");
b1 = load("problem1data\b1.csv");
b2 = load("problem1data\b2.csv");
[m1, n1] = size(A1);
x1 = zeros(n1,1);
[m2, n2] = size(A2);
x2 = zeros(n2,1);
M1 = max(eig(A1'*A1));
M2 = max(eig(A2'*A2));
epsilon = 1e-5;
epsilon_star = 1e-6;
x1_list = [x1];
x2_list = [x2];

flag = 1; 
while flag
    former_x1 = x1;
    temp = x1 - gf(A1, b1, x1)/M1;
    x1(temp>=1/M1) = temp(temp>=1/M1)-1/M1;
    x1(temp<=-1/M1) = temp(temp<=-1/M1)+1/M1;
    x1(-1/M1<temp & temp<1/M1) = 0;
    x1_list = [x1_list x1];
    if(abs(f(A1, b1, x1) - f(A1, b1, former_x1)) < epsilon)
        flag = 0;
    end
end
h1_list = [];
for i=1:size(x1_list,2)
    h1_list = [h1_list h(A1, b1, x1_list(:,i))];
end
while abs(f(A1, b1, x1) - f(A1, b1, former_x1)) > epsilon_star
    former_x1 = x1;
    temp = x1 - gf(A1, b1, x1)/M1;
    x1(temp>=1/M1) = temp(temp>=1/M1)-1/M1;
    x1(temp<=-1/M1) = temp(temp<=-1/M1)+1/M1;
    x1(-1/M1<temp & temp<1/M1) = 0;
end

flag = 1; 
while flag
    former_x2 = x2;
    temp = x2 - gf(A2, b2, x2)/M2;
    x2(temp>=1/M2) = temp(temp>=1/M2)-1/M2;
    x2(temp<=-1/M2) = temp(temp<=-1/M2)+1/M2;
    x2(-1/M2<temp & temp<1/M2) = 0;
    x2_list = [x2_list x2];
    if(abs(f(A2, b2, x2) - f(A2, b2, former_x2)) < epsilon)
        flag = 0;
    end
end
h2_list = [];
for i=1:size(x2_list,2)
    h2_list = [h2_list h(A2, b2, x2_list(:,i))];
end
while abs(f(A2, b2, x2) - f(A2, b2, former_x2)) > epsilon_star
    former_x2 = x2;
    temp = x2 - gf(A2, b2, x2)/M2;
    x2(temp>=1/M2) = temp(temp>=1/M2)-1/M2;
    x2(temp<=-1/M2) = temp(temp<=-1/M2)+1/M2;
    x2(-1/M2<temp & temp<1/M2) = 0;
end

figure(1),
plot(log(1:length(h1_list)), log(sqrt(sum((x1_list - x1).^2)))),
title('1-log(norm(x_k - x^*))-log(k)'),
xlabel('log(k)'),
ylabel('log(norm(x_k - x^*))');
figure(2),
plot(log(1:length(h1_list)), log(h1_list)),
title('1-log(h(x_k))-log(k)'),
xlabel('log(k)'),
ylabel('log(h(x_k))');
figure(3),
plot(log(1:(length(h1_list)-1)), log(h1_list(1:end-1) - h1_list(2:end))),
title('1-log(h(x_{k-1} - h_k))-log(k)'),
xlabel('log(k)'),
ylabel('log(h(x_{k-1} - h_k))');

figure(4),
plot(log(1:length(h2_list)), log(sqrt(sum((x2_list - x2).^2)))),
title('2-log(norm(x_k - x^*))-log(k)'),
xlabel('log(k)'),
ylabel('log(norm(x_k - x^*))');
figure(5),
plot(log(1:length(h2_list)), log(h2_list)),
title('2-log(h(x_k))-log(k)'),
xlabel('log(k)'),
ylabel('log(h(x_k))');
figure(6),
plot(log(1:(length(h2_list)-1)), log(h2_list(1:end-1) - h2_list(2:end))),
title('2-log(h(x_{k-1} - h_k))-log(k)'),
xlabel('log(k)'),
ylabel('log(h(x_{k-1} - h_k))');
%%
%%
function result = h(A, b, x)
    result = norm(A*x - b)^2/2 + norm(x, 1);
end
function result = f(A, b, x)
    result = norm(A*x - b)^2/2;
end
function result = gf(A, b, x)
    result = A'*(A*x-b);
end