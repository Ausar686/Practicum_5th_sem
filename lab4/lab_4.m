%% TASK 1
clc;

A = [1, 1];
B = [0, 0];
C = [-1, 1];
[x1, x2, x3, x4] = biquadsolve(A, B, C);
%[x1, x2] = biquadsolve(A, B, C);
x1
x2
x3
x4
res1 = abs(A.*(x1.^4) + B.*(x1.^2) + C);
res2 = abs(A.*(x2.^4) + B.*(x2.^2) + C);
res3 = abs(A.*(x3.^4) + B.*(x3.^2) + C);
res4 = abs(A.*(x4.^4) + B.*(x4.^2) + C);
res1
res2
res3
res4

clear;
%% TASK 3
clc;

n = 50;
val1 = zeros(1,n);
val2 = zeros(1,n);
val3 = zeros(1,n);
x = linspace(1,n,n);
for i=1:n
    A = rand(i,i);
    [Qm,Rm] = qr_m(A);
    [Q,R] = qr(A);
    [Qc,Rc] = qr_c(A);
    val1(i) = norm(Qm*Rm-A);
    val2(i) = norm(Q*R-A);
    val3(i) = norm(Qc*Rc-A);
end

figure;
plot(x,val1, 'r', x, val2, 'b', x,val3, 'g');
legend('qr_m','qr','qr_c');
xlabel('Размерность матрицы');
ylabel('Величина нормы');
title('Нормы разности');

clear;
%% TASK 4
clc;

n = 50; % количество размерностей
k = 5; % количество "тик/ток"
mTime = zeros(1,n);
cTime = zeros(1,n);
defaultTime = zeros(1,n);
tmp1 = zeros(1,k);
tmp2 = zeros(1,k);
tmp3 = zeros(1,k);
x = linspace(1,n,n);
for i=1:n
    for j=1:k
        A = rand(i,i);
        tic;
        qr_m(A);
        tmp1(i) = toc;
        tic;
        qr(A);
        tmp2(i) = toc;
        tic;
        qr_c(A);
        tmp3(i) = toc;
    end
    mTime(i) = mean(tmp1);
    defaultTime(i) = mean(tmp2);
    cTime(i) = mean(tmp3);
end

plot(x, defaultTime, 'r', x, mTime, 'g', x, cTime, 'b');
xlabel("Размерность матрицы");
ylabel("Время");
title('Среднее время работы в зависимости от размерности');
legend('qr', 'qr_m', 'qr_c');

clear;
%% TASK 5 
clc;

n = 50;
k = 3;
mTime = zeros(1,n);
cTime = zeros(1,n);
defaultTime = zeros(1,n);
tmp1 = zeros(1,k);
tmp2 = zeros(1,k);
tmp3 = zeros(1,k);
x = linspace(1,n,n);
x1 = linspace(1,n,10*n);
for i=1:n
    for j=1:k
        A = rand(i,i);
        tic;
        qr_m(A);
        tmp1(i) = toc;
        tic;
        qr(A);
        tmp2(i) = toc;
        tic;
        qr_c(A);
        tmp3(i) = toc;
    end
    mTime(i) = mean(tmp1);
    defaultTime(i) = mean(tmp2);
    cTime(i) = mean(tmp3);
end

degree = 5;

figure;
p = polyfit(x,mTime,degree); % коэффициенты аппроксимирующего полинома
Y = polyval(p,x1); % значения полинома на сетке х1
plot(x1, Y, x, mTime);
legend('Аппроксимация', 'Функция');
xlabel('Размерность');
ylabel('Время');
title('qr_m');

figure;
p = polyfit(x,defaultTime,degree);
Y = polyval(p,x1);
plot(x1, Y, x, defaultTime);
legend('Аппроксимация', 'Функция');
xlabel('Размерность');
ylabel('Время');
title('qr');

figure;
p = polyfit(x,cTime,degree);
Y = polyval(p,x1);
plot(x1, Y, x, cTime);
legend('Аппроксимация', 'Функция');
xlabel('Размерность');
ylabel('Время');
title('qr_c');

clear;
%% TEST
clc;

clear;
