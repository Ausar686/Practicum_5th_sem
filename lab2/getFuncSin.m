function f = getFuncSin(n, l)
f = @(x) sin(2 .* pi .* n .* x ./ l);
end