function val = func_1_1(x)
% f(x) = cos(x^2 - 4abs(x)))
val = cos(x .* x - 4 * abs(x));
end
