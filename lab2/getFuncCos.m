function f = getFuncCos(n, l)
f = @(x) cos(2 .* pi .* n .* x ./ l);
end