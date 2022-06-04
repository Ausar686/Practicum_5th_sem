clc;
drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0);   


% --------------------------- system 1 ----------------------------------
%   x' = xy - y^3 + x^5, x �� R,
%   y' = x^2 - y^3, y �� R,

v = @(x, y) x.^2 + y.^2;
% v > 0 � ������� U={(x,y) �� R2: |x|, |y| < eps}
% v == 0 �� ������� ������� (����� (0,0))
% dv_dt = -2.*((x-y).^2 + y.^2 < 0 � ������� U
% �� ������� �������� ������� ������� ��������������� ���������
t0 = -10;
tf = 10;
n = 20;
tetta = linspace(0, 2*pi, n);
r = 10;
max1 = 10;

figure;
cla;
axis([-max1 max1 -max1 max1]);
hold on;
[X, Y] = meshgrid(-max1:0.1:max1, -max1:0.1:max1);
Z = v(X, Y);
[~,c] = contourf(X, Y, Z, 20, ':');
c.LineWidth = 0.01;
max_val = max(max(max(Z)));
colorbar

% ------------------------ ������� ������� -------------------------
for i=1:n
    y0 = [r*cos(tetta(i)), r*sin(tetta(i))];
    [t, F] = ode45(@sys1, [t0, tf], y0);
    %plot(F(:, 1), F(:, 2), 'y>-', 'LineSmoothing', 'on');
    col = 0;
    for k=1:2:size(F,1)-1
        if (F(k, 1) < max1 && F(k, 2) < max1)
            col = 5*v(F(k, 1),F(k, 2))/(max_val);
            if col >=1
                col = 1; % ����� �� �������� �� ��������� �����
            end
            a = drawArrow([F(k, 1), F(k+1, 1)], [F(k, 2), F(k+1, 2)]);
            a.Color = [1- col 1 1];
            a.LineWidth = 2;
        end
    end 
end
title('������ �������:x'' = xy - y^3 + x^3, y'' = x^2 - y^3');
legend('���������� �������', '����� ������');
hold off

disp("SYSTEM 1 FINISHED");


% --------------------------- system 2 ----------------------------------
%   x' = x^2+2y^3
%   y' = xy^2
v = @(x, y) (x.^2).*y ;
% v > 0 � ������� U=1-�� ��������
% v == 0 �� ������� ������� 
% dv_dt > 0 � ������� U
% ������� ������� ����������� �� �������� (�.�������)
t0 = -10;
tf = 10;
n = 30;
tetta = linspace(0, 2*pi, n);
r = 20;
max1 = 15;

% ------------------------- ����� ������ -----------------------------
figure;
cla;
axis([-max1 max1 -max1 max1]);
hold on;
[X, Y] = meshgrid(-max1:0.1:max1, -max1:0.01:max1);
Z = v(X, Y);
[~,c] = contourf(X, Y, Z, 20, ':');
c.LineWidth = 0.01;
max_val = max(max(max(Z)));
title('������ �������: x'' = x^2+2y^3, y'' = xy^2');

% ------------------------ ������� ������� -------------------------

for i=1:n
    y0 = [r*cos(tetta(i)), r*sin(tetta(i))];
    [t, F] = ode45(@sys2, [t0, tf], y0);
    col = 0;
    for k=1:2:size(F,1)-1
        if (F(k, 1) < max1 && F(k, 2) < max1)
            col = v(F(k, 1),F(k, 2))/max_val/5;
            if col >=0.5
                col = 0.5;
            end
            if col <=-0.5
                col = -0.5;
            end
            %col
            col = [0.5-col 1 1];
            plot([F(k, 1), F(k+1, 1)], [F(k, 2), F(k+1, 2)], 'Color', col)
            plot (F(k+1, 1), F(k+1, 2), '.', 'Color', col)
        end
    end 
end
legend('���������� �������', '����� ������');
hold off

disp("FINISHED");
clear;

% x' = y-x+xy 
% y' = x-y-x^2-y^3
function dxdt = sys1(t, x)
dxdt = [x(2) - x(1) + x(1).*x(2); ...
    x(1) - x(2) - x(1).^2 - x(2).^3];
end

%   x' = x^2+2y^3
%   y' = xy^2
function dxdt = sys2(t, x)
dxdt = [ x(1).^2 + 2.*x(2)^3; ...
    x(1).*x(2).^2];
end
