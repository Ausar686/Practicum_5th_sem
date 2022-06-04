close all
clc;
n = 2;

f = @(x) cos(x) - x./pi;
f1 = @(x) cos(x);
f2 = @(x) x./pi;
%g = @(x) 0;

x = -10:0.0001:10;
plot(x, f1(x), x, f2(x));
%plot(x, f(x), x, g(x));
hold on
[x1, y1] = ginput(n);
disp(x1);
for i=1:size(x1, 1)
    str = strcat ('-------------- решение № ', num2str(i));
    disp(str)
    answ = fzero (f, x1(i));
    plot(answ, f1(answ), 'b*')
    str = strcat('Корень x = ', num2str(answ));
    str = strcat(str, ", y = ");
    str = strcat(str, num2str(f1(answ)));
    disp(str)
    str = strcat ('Невязка = ', abs(num2str(f2(answ) - f1(answ))));
    disp(str)
end

legend('y = cos(x)', 'y = x/\pi', 'Корни')
hold off;
clear;
