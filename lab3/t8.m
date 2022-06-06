% y'' + y = 2x - pi
% y(0) = 0
% y(pi) = 0
clc;

xmesh = 0:0.1:pi;
analyt = @(x) 2.*x + pi.*cos(x) +1.95*sin(x) - pi;
solinit = bvpinit(xmesh, [-1 2]);
sol = bvp4c(@bvpfcn, @bcfcn, solinit);
figure(2)
plot(sol.x, sol.y(1,:), '-*', xmesh, analyt(xmesh))
legend('численное решение', 'аналитическое решение')

disp('Норма L2')
trapz((abs(sol.y(1, :) - analyt(xmesh))).^2)^0.5
disp('Норма С')
max((abs(sol.y(1, :) - analyt(xmesh))).^2)

clear;

function dydx = bvpfcn(x,y)
dydx = [y(2) 2.*x-pi-y(1)];
end

% Краевые условия
function res = bcfcn(ya,yb)
res = [ya(1) yb(1)];
end
