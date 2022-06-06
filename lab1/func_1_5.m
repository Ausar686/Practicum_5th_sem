function res = func_1_5(x, y, n, m)
% Функция всех пар декартового произведения
x_new = reshape(repmat(x, m, 1), [n * m, 1]); % Делаем столбец из x1.(m раз)..x1x2...x2....xn...xn
y_new = reshape(repmat(y, n, 1)', [n * m, 1]); % Делаем столбец из n штук y1...ymy1y2...ym
res = horzcat(x_new, y_new); % Соединяем их
end
