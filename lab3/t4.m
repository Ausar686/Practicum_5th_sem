% x'' = 0
% x = Cx + C1
% движение внутри окружности
% y = y(x)
%задать кол-во отскоков (event-ов)
clc;
t0 = 0;
x0 = [0.4,0];
v0 = 10 * [3,-1.5];
alpha = 1.001;
num_rebounds = 10;
 
thetta = 0:0.1:2*pi+0.1;
global rad
rad = 10;
x_grid = sin(thetta).*rad;
y_grid = cos(thetta).*rad;

figure(1)
cla
hold on
plot(x_grid, y_grid, 'r')
xlim([-rad-1 rad+1])
ylim([-rad-1 rad+1])


for i = 1:1:num_rebounds
    a = v0 * v0';
    b = v0 * x0';
    c = x0 * x0' - rad^2;
    D = b^2 - a * c;
    t1 = (-b + sqrt(D)) / a;
    x1 = x0 + v0 * t1;
    
    n = x1 / norm(x1);
    l = v0;
    ln = (n * l') * n;
    v1 = (l - 2 * ln) / alpha;
    
    F = @(t, x) v0';
    options = odeset('Events', @bounceEvents); 
    [T, X] = ode45(F, [t0:0.01:t0 + 10000], x0, options);
    for i = 1:size(X,1)-1
        plot(X(i, 1), X(i, 2), 'bo');
        pause(0.01)
        plot(X(i, 1), X(i, 2), 'wo');
    end
    
    x0 = x1;
    v0 = v1;
    t0 = t1;
end
hold off

clear;



function [value,isterminal,direction] = bounceEvents(t,y)
    global rad;
    value = y(1)^2 + y(2)^2 - rad^2;
    isterminal = 1;
    direction = 1;
end
