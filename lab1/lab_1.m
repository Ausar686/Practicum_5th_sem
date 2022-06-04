%% TASK 1
a = input('');
b = input('');
n = input('');
if ~isnumeric(a) || ~isnumeric(b) || ~isnumeric(n) || n <= 1 || mod(n, 1) >= 0.0000001 % �������� ������������ �������
    clear;
    error('Wrong arguments')
end
if a > b % �������������� ������
    t = a;
    a = b;
    b = t;
end
dots = linspace(a, b, n); % ����� �����
values = func_1_1(dots); % �������� � �����
[y_max, max_index] = max(values); % y_max - ��������, max_index - ������ ��������� � �������
[y_min, min_index] = min(values); % ���������� ��� ��������
x_min = dots(min_index); % ������� �������� ����� ���.
x_max = dots(max_index); % ������� �������� ����� ����.
disp("���������� ������� - " + num2str(y_min)); % ���������� ���.(��������)
disp("���������� �������� - " + num2str(y_max)); % ���������� ����.(��������)
plot(dots, values);
xlabel("x");
ylabel("x^2 cos(|x|)");
legend("y = x^2 cos(|x|)");
repmat(y_max, size(x_max)); % ����������� ������� y_max � x_max
repmat(y_min, size(x_min)); % ���������� ��� ���.
text(x_min, y_min, '\leftarrow y_{min}') % �������� �� ������� ����� �������� ����������
text(x_max, y_max, '\leftarrow y_{max}') % ���������� ��� ���������
clear; % ������� ���������� �����
%% TASK 2
n = input('');
disp("�������� �� �������� ����� n ������ ���������: " + num2str(isprime(n)));
disp("����� �����, ��������������� ������� ������...")
if n < 7
    disp("����� ����� ���!");
else
    %x = 1:n / 7:2; % ������ �������� ����� ����� �� 1 �� n/7
    %disp(7 * x); % ������� �����
    disp(7:14:n);
end
res1 = (ones(n) .* (2:(n + 1)))'; % ������� n*n �� ������, ����� �������� �� ������� �� 2,3,4...
disp(res1);
res2 = reshape(1:((n + 1) ^ 2), [n + 1, n + 1])'; % ����� �������� �� 1 �� (n+1)^2, ����� ������ ������ � �������������
disp(res2);
res3 = reshape(res2(n^2 : (n + 1)^2), [n + 1, 2]); % ���� ���������, ����� �������� ����� � ���� ���� �������� � ������� �� � ������ �������
disp(res3);
clear;
%% TASK 3
matrix = round(rand(5, 8) * 21 - 10.5); % ����� ������� � ����. �������. ���������� �� [0, 1],
% ��������� �� 21 (�.�. ��� U[0, 21]) � �������� 10.5 (�.�. U[-10.5, 10.5])
% ����� � ������ ����� ����� ��������� ���������� (n-1/2, n+1/2), �������
% �������, �� ����������� ��������� � ����� ����� 0.
res1 = max(diag(matrix)); % ������������ ������� �� ���������
sums = sum(matrix, 1); % ���������� �����
prods = prod(matrix, 1); % ���������� ������������
res2 = max(sums ./ prods); % ������������ ��������� ���������� ���� � ���������� �������������
res3 = min(sums ./ prods); % ����������� ��������� ���������� ���� � ���������� �������������
new_matrix = sortrows(matrix, 'descend'); % ��������� � �������� ������������������ �������
disp("���������� ������������ ������� - " + num2str(res1));
disp("���������� ��������� ����� � �������������: " + num2str(res2));
disp("���������� ��������� ����� � �������������: " + num2str(res3));
disp(new_matrix);
clear;
%% TASK 4
n = input('');
m = input('');
matrix = rand(n, m);
[r, g, b] = func_1_4(matrix, n, m);
disp("MATRIX = ");
disp(matrix);
disp("R = ");
disp(r);
disp("G = ");
disp(g);
disp("B = ");
disp(b);
clear;
%% TASK 5
n = input('');
m = input('');
x = rand(1, n); %��������� �������� ������ x
y = rand(1, m); %��������� �������� ������ y
res = func_1_5(x, y, n, m);
disp(x);
disp(y);
disp(res);
clear;
%% TASK 6
n = input('');
x = rand(3, n); % ������� ��������� ��������
x_new = reshape(repmat(x, n, 1), [3, n * n]); % ������ ������ ������ 5: ��������� ������������  ��-�� �������� �� ����
y_new = repmat(x, 1, n); 
crosses = cross(x_new, y_new); % �� ��������� ������������ ����� ��������� ������������
abses = reshape(sqrt(sum(crosses .* crosses, 1)), [n, n]); % � �� ��������� ������������ - �� ������. ����� ����� ���� ������� ������� ����������.
disp(abses);
clear;
%% TASK 7
% �������� ������ ������� �� ������� ��������
a1 = [1, 2, 3, 4];
b1 = [-1, -2, -3, -4];
a2 = [1, 2, 3];
b2 = [3, 4, 5];
disp(func_1_7(a1, b1));
disp(func_1_7(a2, b2));
clear;
%% TASK 8
n = input('');
k = input('');
x = rand(k, n); % ������� ��������� ��������
x_new = reshape(repmat(x, n, 1), [k, n * n]); % ������ ������ ������ 5: ��������� ������������  ��-�� �������� �� ����
y_new = repmat(x, 1, n);
diffs = reshape(sqrt(sum((x_new - y_new) .* (x_new - y_new), 1)), [n, n]); % ������ ������ 6 ��� ����������
disp(diffs);
clear;
%% TASK 9
% ��� ���������� ������ ���� ����� ������������� ������� ������� ������
% ������� ������ ��� ���������� ������
n_times = 100;
max_dim = 300;
slow_average = zeros(1, max_dim); % ������� ����� ������ ����� ������� �� �������� ���������� �������
fast_average = zeros(1, max_dim); % ������� ����� ������ ����������� ������ �� �������� ���������� �������
for n = 1:max_dim
    slow_res = zeros(1, n_times); % ���������� ������ ����� ������� �� �������� ������� �������
    fast_res = zeros(1, n_times); % ���������� ������ ����������� ������ �� �������� ������� �������
    for i = 1:n_times
        a = rand(n);
        b = rand(n);
        tic();
        func_1_9(a, b, n, n);
        slow_res(i) = toc();
        tic();
        a + b;
        fast_res(i) = toc();
    end
    slow_average(n) = mean(slow_res);
    fast_average(n) = mean(fast_res);
end
plot(1:max_dim, slow_average, 1:max_dim, fast_average);
xlabel("����������� �������");
ylabel("������� ����� ���������� (�)");
legend("������������ ��������", "���������� �����");
clear;
%% TASK 10
% �������� ������ ������� �� ������� �������
a1 = [1, 2, 3, 2, 1];
a2 = [1, 2, 3, 4, 1];
a3 = (1);
a4 = [2, 2];
a5 = [6, 9];
disp(func_1_10(a1)); % 1
disp(func_1_10(a2)); % 0
disp(func_1_10(a3)); % 1
disp(func_1_10(a4)); % 1
disp(func_1_10(a5)); % 0
clear;
%% TASK 11
a = input('');
if a <= 0
    clear;
    error("�������� ������� �������");
end
b = input('');
n = input('');
x = rand(1, n) * a; % n ��������� ������� ���������� �������������� �� [0, a]
num = sum((x >= b)) / n; % ������� ���� �����, ������� b
if num <= a / (2 * b)
    disp("����������� � ���������!");
else
    disp("���-�� ����� �� ���...");
end
%disp(x);
disp(num);
disp(a / (2 * b));
clear;
%% TASK 12
% ����� ������� ����� ������������� exp(-x^2), ������� � ����� ����� ����� 0: F(a)= 0
% ������� ����������
a = -2; % ������� ����� ����� �������
b = 2; % ������� ������ ����� �������
n = 30; % ������� ����� �������� ���������
m = n; % ���������� �������� ��������� ���������� [x_k, x_k+1]
h = (b - a) / n; % ������ �������� ����
gauss = @(x) exp(-x .* x); % ��������� ������� exp(-x^2)
tiledlayout(2, 2);
%������� ���������������
[x_rect, y_rect] = rectangles(a, b, n, m, gauss);
% ������� ��������
[x_sim, y_sim] = simpson(a, b, n, m, gauss);
% ������� ��������
x_trap = linspace(a, b, n + 1);
y_trap = zeros(1, n);
for k = 1:n
    x_k = linspace(a, a + k * h, m * k + 1);
    y_k = gauss(x_k);
    y_trap(k) = trapz(x_k, y_k);
end;
y_trap = horzcat(0, y_trap);
% ����� ������
ax1 = nexttile;
plot(ax1, x_rect, y_rect, x_sim, y_sim, x_trap, y_trap);
xlabel("x");
ylabel("������������� exp(-x^2)");
legend("����� ���������������", "����� ��������", "����� ��������");
% ����� ������ �������
ax2 = nexttile;
n_times = 100;
n_dots = 100;
rect_times = zeros(1, n_dots);
sim_times = zeros(1, n_dots);
trap_times = zeros(1, n_dots);
rect_cur = zeros(1, n_times);
sim_cur = zeros(1, n_times);
trap_cur = zeros(1, n_times);
for k = 1:n_dots
    for i = 1:n_times
        tic();
        rectangles(a, b, n, k, gauss);
        rect_cur(i) = toc();
        tic()
        simpson(a, b, n, k, gauss);
        sim_cur(i) = toc();
        tic()
        x_k = linspace(a, a + h, k + 1);
        y_k = gauss(x_k);
        trapz(x_k, y_k);
        trapz_cur(i) = (n + 1) * toc();
    end
    rect_times(k) = mean(rect_cur);
    sim_times(k) = mean(sim_cur);
    trap_times(k) = mean(trap_cur);
end
plot(ax2, 1:n_dots, rect_times, 1:n_dots, sim_times, 1:n_dots, trap_times);
xlabel("���������� �������� ���������");
ylabel("������� ����� �� ���������� (�)");
legend("����� ���������������", "����� ��������", "����� ��������");
% ��������� �������� ���������� ����������
ax3 = nexttile;
n_dots = 100;
n_start = 3;
rect_diff = zeros(1, n_dots);
sim_diff = zeros(1, n_dots);
trap_diff = zeros(1, n_dots);
for k = n_start:n_dots
    x_tmp1 = linspace(a, b, k + 1);
    y_tmp1 = gauss(x_tmp1);
    x_tmp2 = linspace(a, b, 2 * k + 1);
    y_tmp2 = gauss(x_tmp2);
    [tmp, rect_k] = rectangles(a, b, k, 1, gauss);
    [tmp, rect_2k] = rectangles(a, b, 2 * k, 1, gauss);
    [tmp, sim_k] = simpson(a, b, k, 1, gauss);
    [tmp, sim_2k] = simpson(a, b, 2 * k, 1, gauss);
    trap_k = trapz(x_tmp1, y_tmp1);
    trap_2k = trapz(x_tmp2, y_tmp2);
    rect_diff(k) = abs(rect_k(end) - rect_2k(end));
    sim_diff(k) = abs(sim_k(end) - sim_2k(end));
    trap_diff(k) = abs(trap_k - trap_2k);
end;
rect_diff = rect_diff(n_start:n_dots);
sim_diff = sim_diff(n_start:n_dots);
trap_diff = trap_diff(n_start:n_dots);
plot(ax3, n_start:n_dots, rect_diff, n_start:n_dots, sim_diff, n_start:n_dots, trap_diff);
xlabel("���������� �������� ���������");
ylabel("������ �������� �� n � 2n ��������");
legend("����� ���������������", "����� ��������", "����� ��������");
clear;
%% TASK 13
n = 100; % ����� ��������� ��������� ����
x = ones(1, n); % x = 1 ��� ������� �� n ��������� ����
h = logspace(-10, -1, n); % �������� ����
abs_right = abs(func_1_13_2(x) - (func_1_13_1(x + h) - func_1_13_1(x)) ./ h); % ������ ������� �������� � ������ ���������� �����������
abs_center = abs(func_1_13_2(x) - (func_1_13_1(x + h) - func_1_13_1(x - h)) ./(2 * h)); % ������ ������� �������� � ����������� ���������� �����������
loglog(h, abs_center, h, abs_right);
xlabel("������ ����");
ylabel("������ ���������� �� �����������");
legend("�����. ����. ��-��", "����. ����. ��-��");
clear;
%% TEST
square = @(x) x .* x;

clear;