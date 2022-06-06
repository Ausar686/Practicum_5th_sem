close all
clc;
a = 10;

a1 = 0.32; % 1/pi
nAnswers = 101;
x = linspace(-a1, a1, nAnswers*2);
xAns = zeros(1, size(x, 1));
for i=1:nAnswers*2
    xAns(i) = fzero(@myf, x(i));
end
figure(1)
x(1) = -a;
x(end) = a;
plot(x, xAns)
title('Ближайший корень функции f(x) = 0')

%График функции может быть не похож на то, что есть на самом деле из-за
%малости n: мы всегда попадаем в ту часть, где f < 0 (в окрестности x = 0)
figure(2)
hold on
nDots = 10001;
x = linspace(-a, a, nDots);
plot(x, myf(x), [-a a], [0 0])
plot(xAns, zeros(size(xAns, 2)), '*')
legend('Функция', 'Оси', 'Корни')
hold off

clear;

function y = myf (x)
    y = x.*sin(1./(x));
    if (length(x) == 1)&&(x == 0)
        y = 0;
    elseif find(x == 0)
        y(find(x == 0)) = 0;
    end
end
