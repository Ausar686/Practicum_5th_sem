function convergenceFunc(fn, f, a, b, n, convType)
fig = figure;
nDots = 101;
x = linspace(a, b, nDots);
y = f(x);
mov(1:n) = struct('cdata', [], 'colormap', []);
for i = 1:n
    pause(0.2);
    plot(x, fn(i, x), x, y);
    legend(func2str(fn), func2str(f));
    xlabel('x');
    ylabel('y');
    if strcmp(convType, 'Pointwise')
        title('');
    end
    if strcmp(convType, 'Uniform')
        title(string(max(abs(y - fn(i, x)))));
    end
    if strcmp(convType, 'L2')
        title(string(sqrt(trapz(x, (y - fn(i, x)).^2))));
    end
    mov(i) = getframe(fig);
end
%movie(mov);
end