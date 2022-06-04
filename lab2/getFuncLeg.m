function f = getFuncLeg(n)
if n == 0
    f = @(x) 1;
elseif n == 1
    f = @(x) x;
else
    P1 = getFuncLeg(n - 1);
    P2 = getFuncLeg(n - 2);
    f = @(x) (2 * n - 1) .* x .* P1(x) ./ n - (n - 1) .* P2(x) ./ n;
end
end