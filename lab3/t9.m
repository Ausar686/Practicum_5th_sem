% минимум ф-ии многих переменных методом Ньютона
% ф-ия, градиент, гессиан, начальное приближение

% для 2х переменных построить набор линий уровня, на них отметить шаги
% алгоритма. Сравнить с fminbd
clc;

f = @(x, y) x.^2+ 5*y.^2 - sin(x.^2+y) + exp(x);
grad = @(x, y) [2.*x - 2.*x.*cos(x.^2+y) + exp(x); 10.*y - cos(x.^2+y)];
x1 = [1; 1];
x0 = x1;
eps = 0.001;

x = -5:0.1:5;
y = -5:0.1:5;
[x, y] = meshgrid(x, y);
figure(1)
subplot(1,2, 1)
hold on
contour(x, y, f(x, y),30)
xlim([-5 5])
ylim([-5 5])
a = 0.1
title('x^2+ 5*y^2 - sin(x^2+y) + exp(x)')
while sum(grad(x0(1), x0(2)).^2) > eps^2
    x0 = x1;
    plot (x1(1), x1(2), 'r*')
    x1 = x0 - a*grad(x0(1), x0(2));
    %pause(1);
    sum(grad(x0(1), x0(2)).^2)
    pause(0.1)
end

plot (x1(1), x1(2), 'r*')
min_val = inf;
min_x = inf;
min_y = inf;
for i=-5:0.1:5
    f1 = @(x) x.^2+ 5*i.^2 - sin(x.^2+i) + exp(x);
    xmin = fminbnd(f1, -1, 1)
    val = f(xmin, i);
    if val < min_val
        min_val = val;
        min_x = xmin;
        min_y = i;
    end
end

%clear;
%% функция 2х переменных
clc;

f = @(x) xmin.^2+ 5*x.^2 - sin(xmin^2+x);
ymin = fminbnd(f, -1, 1)
plot (xmin, ymin, 'o')
hold off


f = @(x, y) -1./(x.^2+1) + sin(y) + 1;
grad = @(x, y) [2.*x./(x.^2+1).^2; cos(y)];
x1 = [0; 4];
x0 = x1;
eps = 0.001;

x = -2:0.1:2;
y = 2:0.1:6;
[x, y] = meshgrid(x, y);
subplot(1,2,2)
hold on
contour(x, y, f(x, y),30)
xlim([-2 2])
ylim([2 6])
title('-1/(x^2+1) + sin(y) + 1')
while sum(grad(x0(1), x0(2)).^2) > eps^2
    x0 = x1;
    plot (x1(1), x1(2), 'r*')
    x1 = x0 - a*grad(x0(1), x0(2));
    f(x1(1), x1(2))
    pause(0.1);
end

plot (x1(1), x1(2), 'r*')
min_val = inf;
min_x = inf;
min_y = inf;
for i=-5:0.1:5
    f1 = @(x) -1./(x.^2+1) + sin(i) + 1;
    xmin = fminbnd(f1, -1, 1)
    val = f(xmin, y);
    if val < min_val
        min_val = inf;
        min_x = inf;
        min_y = inf;
    end
end

plot (xmin, ymin, 'go')
hold off

%clear;

%% функция 3х переменных

f3 = @(x) x(1).^2 + 2.*x(2)^2 + x(1)^2*x(2)^2 +x(3)^2 + exp(x(2).^2 +x(3).^2) - x(2) + x(3);
grad3 = @(x, y, z) [2*x + 2*x.*y^2;  4*y + 2*x.^2*y + 2*y.*exp(y.^2 +z.^2) - 1;  2*z + 2*z.*exp(y.^2 +z.^2) + 1];
x1 = [1; 1; 1];
x0 = [1; 1; 1];
eps = 0.001;
n = 100;

while sum(grad3(x0(1), x0(2), x0(3)).^2) > eps.^2
    x0 = x1;
    x1 = x0 - a*grad3(x0(1), x0(2), x0(3))
end
x1
z = fminsearch(f3, [1,1,1])

min_val = inf;
min_x = inf;
min_y = inf;
min_z = inf;
for i=-1:0.01:1
    for j=-1:0.01:1
        f1 = @(x) x.^2 + 2.*i^2 + x.^2*i^2 +j^2 + exp(i.^2 +j.^2) - i + j;
        xmin = fminbnd(f1, -1, 1);
        val = f3([xmin i j]);
        if val < min_val
            min_val = val;
            min_x = xmin;
            min_y = i;
            min_z = j;
        end
    end
end
min_x
min_y
min_z
%pause
%% Для 1 переменной: сравнение с fminbnd
clc;

f = @(x) sin(x);
grad = @(x) cos(x);
a = 1;
x = 1; % начальное приближение
while (abs(grad(x)) > 0.01)
    x = x - a .* grad(x);
    a = a ./ (a + 1);
end
x1 = -pi;
x2 = pi;
ans = fminbnd(f, x1, x2);
disp(abs(x - ans));
x
ans
clear;