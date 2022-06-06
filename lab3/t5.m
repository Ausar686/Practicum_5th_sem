% dx = alpha * x - gamma * xy
% dy = -beta * y + delta * xy
%close all
clc;
global alpha beta gamma delta
alpha = 0.1; % коэф рождаемости жертв
beta = 2; % коэф смертности жертв за счет поедания хищниками
gamma = 0; % коэф смертности хищников
delta = 0; % коэф рождаемости хищников за счёт поедания жертв

x0 = 1000; % начальная популяция жертв
y0 = 20;  % начальная популяция хищников
t0 = 0;
t1 = 5;
tspan = t0:0.01:t1;

analytic_x = @(t) x0.*exp(alpha.*t); % аналитическое решение
analytic_y = @(t) y0.*exp(-beta.*t); 
[t, y] = ode45(@f1, [t0, t1], [x0; y0]);

figure;
n = size(y, 1);
hold on
plot(y(:, 1), y(:, 2), 'y');
plot(analytic_x(tspan), analytic_y(tspan), 'g');
xlabel('размер популяции жертв')
ylabel('размер популяции хищников')
title('фазовые кривые')

grid on
n_dots = min(n, 3000);
for i = 1:n_dots
    plot(y(i, 1), y(i, 2), 'ro');
    pause(0.0001);
end
legend('численное решение', 'аналитическое решение')
hold off

figure;
plot3(y(:, 1), y(:, 2), t, 'y');
xlabel('размер популяции жертв')
ylabel('размер популяции хищников')
zlabel('время t')
title('интегральные кривые')
hold on
plot3(analytic_x(tspan), analytic_y(tspan), tspan, 'g');
grid on
n_dots = min(size(t, 1), 3000);
for i = 1:n_dots
    plot3(y(i, 1), y(i, 2), t(i), 'ro');
    pause(0.0001);
end
legend('численное решение', 'аналитическое решение')
hold off
clear;

    
 function dxdt = f1(~, x)
    global alpha beta gamma delta
    x1 = x(1);
    y1 = x(2);
    dxdt = [alpha.*x1 - gamma.*x1.*y1; -beta.*y1 + delta.*x1.*y1];
 end
