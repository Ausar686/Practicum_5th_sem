function viewPossible(points, V, L)
figure;
lbx = -40;
ubx = 40;
xStep = 1;
lby = -40;
uby = 40;
yStep = 1;
x = lbx:xStep:ubx;
y = lby:yStep:uby;
N = length(x);
[X, Y] = meshgrid(x, y);
Z = zeros(N);
dots = ones(N, 2);
dots(:, 2) = y;
for i = 1:N
    dots(:, 1) = x(i);
    d = pdist2(dots, points);
    s = V ./ (1 + d);
    signal = sum(s, 2);
    Z(:, i) = signal;
end
m = contourf(X, Y, Z, [L, L], 'black-');
hold on;
plot(points(:, 1), points(:, 2), 'red*');
xlabel('x');
ylabel('y');
hold off;
if m(2) >= (size(m, 2) - 1)
    disp("Односвязное");
else
    disp("Не односвязное / область выходит за границы");
end