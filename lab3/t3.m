clc;
t0 = 0; t1 = 1;

A = @(t) [1 1; 1 1];
%Через expm
disp('1. expm')
expm(A(t1))

%Через ряд
Aexp = @(t) [1 0; 0 1];
Atmp = @(x) [1 0; 0 1];
for i=1:30
    Atmp = @(x) Atmp(x)*A(x)/i;
    Aexp = @(x) Aexp(x) + Atmp(x);
end
disp('2. Сумма ряда')
Aexp(t1)


% Численное решение задачи коши для ОДУ
% Y'(t) = A*Y(t)
% Y(0) = E
disp('3. ode45')
opts = odeset('Vectorized','on');
[t, y] = ode45(@(t, y) reshape(A(t1)*reshape(y, 2, 2), 4, 1), [t0,t1], [1, 0, 0, 1], opts);
reshape(y(end, :), 2, 2)
clear,