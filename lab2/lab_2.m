%% TASK 1
clc;
f = @(x) x.* x;
grid = linspace(0, 1, 11);
smallGrid = linspace(0, 1, 21);
compareInterp(grid, smallGrid, f);
clear;
%% TASK 2
clc;
antiNearestFunction = @(x) x .^ 2;
antiLinearFunction = @(x) sin(200 .* pi .* x);
antiSplineFunction = @(x) sin(20 .* pi .* x);
antiPchipFunction = @(x) sin(10 .* pi .* x);
nearestFunction = @(x) round(x);
linearFunction = @(x) x .^ 3;

gridNearest = linspace(0, 1, 11);
smallGridNearest = linspace(0, 1, 201);
gridLinear = linspace(0, 1, 11);
smallGridLinear = linspace(0, 1, 1001);
gridSpline = linspace(0, 1, 11);
smallGridSpline = linspace(0, 1, 201);
gridPchip = linspace(0, 1, 11);
smallGridPchip = linspace(0, 1, 201);
goodGridNearest = linspace(0, 10, 101);
goodSmallGridNearest = linspace(0, 10, 1001);

valuesBadNearest = interp1(gridNearest, antiNearestFunction(gridNearest), smallGridNearest, 'nearest');
valuesBadLinear = interp1(gridLinear, antiLinearFunction(gridLinear), smallGridLinear, 'linear');
valuesBadSpline = interp1(gridSpline, antiSplineFunction(gridPchip), smallGridSpline, 'spline');
valuesBadPchip = interp1(gridPchip, antiPchipFunction(gridPchip), smallGridPchip, 'pchip');
valuesGoodNearest = interp1(goodGridNearest, nearestFunction(goodGridNearest), goodSmallGridNearest, 'nearest');
valuesGoodLinear = interp1(gridLinear, linearFunction(gridLinear), smallGridLinear, 'linear');

tiledlayout(3,2);
tile1 = nexttile;
plot(tile1, smallGridNearest, antiNearestFunction(smallGridNearest), smallGridNearest, valuesBadNearest);
legend("x^2", "Nearest-интерполяция");
tile2 = nexttile;
plot(tile2, smallGridLinear, antiLinearFunction(smallGridLinear), smallGridLinear, valuesBadLinear);
legend("sin(200pi * x)", "Linear-интерполяция");
tile3 = nexttile;
plot(tile3, smallGridSpline, antiSplineFunction(smallGridSpline), smallGridSpline, valuesBadSpline);
legend("sin(20pi * x)", "Spline-интерполяция");
tile4 = nexttile;
plot(tile4, smallGridPchip, antiPchipFunction(smallGridPchip), smallGridPchip, valuesBadPchip);
legend("sin(10pi * x)", "Pchip-интерполяция");
tile5 = nexttile;
plot(tile5, goodSmallGridNearest, nearestFunction(goodSmallGridNearest), goodSmallGridNearest, valuesGoodNearest);
legend("round(x)", "Nearest-интерполяция");
tile6 = nexttile;
plot(tile6, smallGridLinear, linearFunction(smallGridLinear), smallGridLinear, valuesGoodLinear);
legend("x^3", "Linear-интерполяция");
clear;
%% TASK 3
clc;
functionGood = @(x) sin(10 .* pi .* x);
functionBad = @(x) sin(100 .* pi .* x);
maxAbsDerivative2Good = 100 .* (pi .^ 2);
maxAbsDerivative2Bad = 10000 .* (pi .^ 2);
grid = linspace(0, 1, 11);
smallGrid = linspace(0, 1, 101);
absDifferenceGood = abs(functionGood(smallGrid) - interp1(grid, functionGood(grid), smallGrid));
absDifferenceBad = abs(functionBad(smallGrid) - interp1(grid, functionBad(grid), smallGrid)); % Тут везде нули
maxSupposedDifferenceGood = repmat(maxAbsDerivative2Good .* (grid(2) .^ 2) ./ 8, [1, size(smallGrid, 2)]);
maxSupposedDifferenceBad = repmat(maxAbsDerivative2Bad .* (grid(2) .^ 2) ./ 8, [1, size(smallGrid, 2)]); % Тут тоже везде нули
tiledlayout(1, 2);
tile1 = nexttile;
plot(tile1, smallGrid, absDifferenceGood, smallGrid, maxSupposedDifferenceGood);
title("Модуль отклонения для y = sin(10pi * x)");
xlabel("Значение переменной");
ylabel("Абсолютная величина отклонения");
legend("Реальный результат", "Априорный результат");
tile2 = nexttile;
plot(tile2, smallGrid, absDifferenceBad, smallGrid, maxSupposedDifferenceBad);
title("Модуль отклонения для y = sin(100pi * x)");
xlabel("Значение переменной");
ylabel("Абсолютная величина отклонения");
legend("Реальный результат", "Априорный результат");
clear;
%% TASK 4
clc;
fn = @(n, x) x .^n;
f = @(x) 0;
a1 = 0;
b = 1;
n = 20;
convType = 'Pointwise';
convergenceFunc(fn, f, a1, b, n, convType);
clear;
%% TASK 5
clc;
f = @(x) x;
a1 = 0;
b = 10;
n = 20;
meth = 'Trigonometry';
fourierApprox(f, a1, b, n, meth);
clear;
%% TASK 6
clc;
%PART 1
a1 = 0;
b = 1;
nDots = 10001;
f = @(x) sin (10 .* pi .* x);
x = linspace(a1, b, nDots);
y = f(x);

%PART 2
yMax = max(y);
maxIdx = find(y == yMax);
xMax = x(maxIdx(1));
minsY = zeros(1, 0);
minsX = zeros(1, 0);
idxMin = 1;
if y(1) <= y(2)
    minsY(1) = y(1);
    minsX(idxMin) = x(1);
    idxMin = idxMin + 1;
end
for i = 2:nDots - 1
    if (y(i - 1) >= y(i)) && (y(i) <= y(i + 1))
        minsY(idxMin) = y(i);
        minsX(idxMin) = x(i);
        idxMin = idxMin + 1;
    end
end
if y(nDots - 1) >= y(nDots)
    minsY(idxMin) = y(nDots);
    minsX(idxMin) = x(nDots);
    idxMin = idxMin + 1;
end
plot(x, y, minsX, minsY, 'r*', xMax, yMax, 'g*');
hold on;
minIdxLeft = find(minsX < xMax);
minIdxRight = find(minsX > xMax);
distLeft = -1;
distRight = -1;
if ~isempty(minIdxLeft)
    distLeft = min(xMax - minsX(minIdxLeft));
end
if ~isempty(minIdxRight)
    distRight = min(minsX(minIdxRight) - xMax);
end
if (distLeft < 0) && (distRight < 0)
    %Никуда не летим
    disp('NO MINIMUM SPOTTED');
elseif (distLeft < 0) || (distLeft >= distRight && distRight > 0)
    %Летим вправо
    minXRight = minsX(minIdxRight);
    xEndIdx = find(minXRight == xMax + distRight);
    xEnd = minXRight(xEndIdx);
    yEnd = f(xEnd);
    x_Idx = find(x >= xMax);
    x_ = x(x_Idx);
    xWayIdx = find(x_ <= xEnd);
    xWay = x_(xWayIdx);
    yWay = f(xWay);
    comet(xWay, yWay);
    hold off;
elseif (distRight < 0) || (distRight > distLeft && distLeft > 0)
    %Летим влево
    minXLeft = minsX(minIdxLeft);
    xEndIdx = find(minXLeft == xMax - distLeft);
    xEnd = minXLeft(xEndIdx);
    yEnd = f(xEnd);
    x_Idx = find(x <= xMax);
    x_ = x(x_Idx);
    xWayIdx = find(x_ >= xEnd);
    xWay = x_(xWayIdx);
    yWay = f(xWay);
    comet(flip(xWay), flip(yWay));
    hold off;
end
clear;
%% TASK 7
clc;
getEqual(@(t) sin(t), @(t) cos(t), 0, 2 * pi, 10);
getEqual(@(t) sin(3*t), @(t) cos(2*t), 0, 5, 10);
getEqual(@(t) sin(t), @(t) cos(2*t), 0, 10, 30);
clear;
%% TASK 8
clc;
%Эллипс с центром в (0, 0) и полуосями 5 и 4
x0 = [0, 0];
a1 = 5;
a2 = 4;
f = @(x) ((x(1) - x0(1)) ./ a1) .^ 2 + ((x(2) - x0(2)) / a2) ^ 2 - 1;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
N = 15;
drawSet(supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub), N);
%Квадрат с центром в (0, 0) и сторонами, равными 2, и параллельными Oxy
g = @(x) 0;
lb = [-1, -1];
ub = [1, 1];
N = 100;
drawSet(supportLebesgue(g, x0, A, b, Aeq, beq, lb, ub), N);
%Ромб с центром в (0, 0) и полуосями 4 и 3, параллельными Oxy
a1 = 4;
a2 = 3;
nDots = @(x) abs(x(1)) / a1 + abs(x(2)) / a2 - 1;
lb = [];
ub = [];
N = 15;
drawSet(supportLebesgue(nDots, x0, A, b, Aeq, beq, lb, ub), N);
clear;
%% TASK 9
clc;
f = @(x) x(1)^2 + x(2)^2 - 1; % f(x) <= 0 - единичный круг с центром в нуле
x0 = [0, 0]; % начальная точка в центре круга, точке (0, 0)
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
supportFunction = supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub);
vector = [2, 2];
res = supportFunction(vector);
disp(res);
clear;
%% TASK 10
clc;

params.a = 5;
params.b = 4;
params.t = [3, 5];
params.c = [0, 0];
N = 100;

drawPolar(@rhoEllipse, N, params);
drawPolar(@rhoRhombus, N, params);

clear;
%% TASK 11, TASK 12
clc;
figure;
a = 1;
b = 10;
n = 200;
nDots = 31;
nLevels = 20;
x = linspace(a, b, nDots);
y = linspace(a, b, nDots);
[X, Y] = meshgrid(x, y);
Z = sin(X) + cos(Y);
C = Z;
contour(X, Y, Z, nLevels);
title('Линии уровня');
colorbar;

fig = figure;
surface = surf(X,Y,Z,C);
colorbar;
hold on;

zMax = Z(2:nDots-1, 2:nDots-1) > Z(2:nDots-1, 1:nDots-2); % точка слева 
zMax = zMax & (Z(2:nDots-1, 2:nDots-1) > Z(1:nDots-2, 2:nDots-1)); % точка сверху
zMax = zMax & (Z(2:nDots-1, 2:nDots-1) > Z(2:nDots-1, 3:nDots)); % точка справа
zMax = zMax & (Z(2:nDots-1, 2:nDots-1) > Z(3:nDots, 2:nDots-1)); % точка снизу
zMin = Z(2:nDots-1, 2:nDots-1) < Z(2:nDots-1, 1:nDots-2); % точка слева
zMin = zMin & (Z(2:nDots-1, 2:nDots-1) < Z(1:nDots-2, 2:nDots-1)); % точка сверху
zMin = zMin & (Z(2:nDots-1, 2:nDots-1) < Z(2:nDots-1, 3:nDots)); % точка справа
zMin = zMin & (Z(2:nDots-1, 2:nDots-1) < Z(3:nDots, 2:nDots-1)); % точка снизу
[iMin, jMin] = find(zMin);
[iMax, jMax] = find(zMax);
iMax = iMax + 1;
jMax = jMax + 1;
iMin = iMin + 1;
jMin = jMin + 1;
pointsMax = plot3(diag(X(iMax, jMax)), diag(Y(iMax, jMax)), diag(Z(iMax, jMax)), 'red *');
pointsMin = plot3(diag(X(iMin, jMin)), diag(Y(iMin, jMin)), diag(Z(iMin, jMin)), 'green *');

j = 0;
delta = 0.04;
for i = 1:n
    j = j + 1;
    M(j) = getframe(fig);
    surface.ZData = sin(X + delta) + cos(Y + delta);
    surface.CData = surface.ZData;
    %Ищем индексы максимумов
    zMax = surface.ZData(2:nDots-1, 2:nDots-1) > surface.ZData(2:nDots-1, 1:nDots-2); % точка слева
    zMax = zMax & (surface.ZData(2:nDots-1, 2:nDots-1) > surface.ZData(1:nDots-2, 2:nDots-1)); %  точка сверху
    zMax = zMax & (surface.ZData(2:nDots-1, 2:nDots-1) > surface.ZData(2:nDots-1, 3:nDots)); % точка справа
    zMax = zMax & (surface.ZData(2:nDots-1, 2:nDots-1) > surface.ZData(3:nDots, 2:nDots-1)); % точка снизу
    [iMax, jMax] = find(zMax);
    %Ищем индексы минимумов
    zMin = surface.ZData(2:nDots-1, 2:nDots-1) < surface.ZData(2:nDots-1, 1:nDots-2); % точка слева
    zMin = zMin & (surface.ZData(2:nDots-1, 2:nDots-1) < surface.ZData(1:nDots-2, 2:nDots-1)); % точка сверху
    zMin = zMin & (surface.ZData(2:nDots-1, 2:nDots-1) < surface.ZData(2:nDots-1, 3:nDots)); % точка справа
    zMin = zMin & (surface.ZData(2:nDots-1, 2:nDots-1) < surface.ZData(3:nDots, 2:nDots-1)); % точка снизу
    [iMin, jMin] = find(zMin);
    % Начинали с отступом, так что нужно вернуться к старым индексам
    iMax = iMax + 1; 
    jMax = jMax + 1;
    iMin = iMin + 1;
    jMin = jMin + 1;
    %Находим координаты супремумов
    pointsMax.XData = diag(surface.XData(iMax, jMax));
    pointsMax.YData = diag(surface.YData(iMax, jMax));
    pointsMax.ZData = diag(surface.ZData(iMax, jMax)); 
    pointsMin.XData = diag(surface.XData(iMin, jMin));
    pointsMin.YData = diag(surface.YData(iMin, jMin));
    pointsMin.ZData = diag(surface.ZData(iMin, jMin));
    
    drawnow;
    delta = delta + 0.04;
    title('Анимация z = cos(x + \deltat) + sin(y - \deltat)');
end
video = VideoWriter('Animation');
video.FrameRate = 60;
open(video);
for j = 1:length(M)
    writeVideo(video, M(j));
end
close(video);

clear;
%% TASK 13
clc;

points = [-10, 5; 20, -2; 5, -10; -5, 7];
viewPossible(points, 1, 0.2);

points = [-30, 10; 30, -10; 0, 0];
viewPossible(points, 1, 0.15);

clear;
%% TASK 14
clc;

params.ax = -5;
params.bx = 5; 
params.ay = -5;
params.by = 5;
params.az = -5;
params.bz = 5;
params.nDots = 100;
params.faceColor = 'blue';
params.edgeColor = 'none';
alpha = 2;
level = 1;

drawBall(alpha, level, params)

clear;
%% TASK 15
clc;

alphas = [4, 2, 1, 0.5, Inf, 2, 4, 1, 5, 6];
colors = ["red", "green", "yellow", "magenta", "blue"];
edges = ["black", "None", "None", "None", "None"];
drawManyBalls(alphas, colors, edges);

clear;
%% TEST
clc;

clear;
%% TEST 2
clc;

clear;
