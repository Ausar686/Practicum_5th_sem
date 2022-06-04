function [res] = rhoRhombus(params) 
% c - точка центра
a = params.a/2;
b = params.b/2;
c = params.c; 
t = params.t;

if (a <= 0) || (b <= 0)
    error('(a <= 0) or (b <= 0)');
end
val =  max(abs(t(1)) * b, abs(t(2)) * a) + t(1) * c(1) + t(2) * c(2);
point = [c(1) + (b*abs(t(1))>=a*abs(t(2)))*b*sign(t(1)), c(2) + (b*abs(t(1))<a*abs(t(2)))*a*sign(t(2))];

res = [point, val];
end
