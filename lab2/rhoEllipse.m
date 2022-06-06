function [res] = rhoEllipse(params) 
% c - точка центра
% alpha -  угол поворота
a = params.a;
b = params.b;
t = params.t;
c = params.c;

if (a <= 0) || (b <= 0)
    error('(a <= 0) || (b <= 0)');
end
T = [a*a, 0; 0, b*b];
Tt = T*t';
point = [c(1) c(2)] + Tt'./sqrt(t(1)*Tt(1) + t(2)*Tt(2)); % значение опорной точки
val = sqrt(a^2*t(1)^2 + b^2*t(2)^2) + dot(t, [c(1) c(2)]); %значение опорной функции

res = [point, val];
end
