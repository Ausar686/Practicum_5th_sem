function res = func_1_9(a, b, n, m)
res = zeros(n, m);
for i = 1:n
    for j = 1:m
        res(i, j) = a(i, j) + b(i, j);
    end
end
end