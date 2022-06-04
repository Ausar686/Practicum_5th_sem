function fourierApprox(f, a, b, n, meth)
nDots = 101;
x = linspace(a, b, nDots);
y = f(x);
len = b - a;
y_i = zeros(1, length(y));
movie(1:n) = struct('cdata', [], 'colormap', []);
if strcmp(meth, 'Trigonometry')
    a_0 = (2 ./ len) .* trapz(x, y);
    y_i = y_i + a_0 ./ 2;
    for i = 1:n
        pause(0.2);
        f_i = getFuncCos(i, len);
        g_i = getFuncSin(i, len);
        a_i = (2 ./ len) .* trapz(x, y .* f_i(x));
        b_i = (2 ./ len) .* trapz(x, y .* g_i(x));
        y_i = y_i + a_i .* f_i(x) + b_i .* g_i(x);
        plot(x, y, x, y_i);
        legend('f(x)', 'i-ая частичная сумма ряда');
        xlabel('x');
        ylabel('y');
        title('Тригонометрическая система функций');
        mov(i) = getframe();
    end
end
if strcmp(meth, 'Chebyshev')        
    x_tmp = linspace(-1+0.0001, 1-0.0001, nDots);
    y_i = 1 ./ pi * trapz(x_tmp, y ./ (1 - x_tmp .^ 2) .^ (0.5));
    j = 0;
    for i = 1:n
        f_i = getFuncCheb(i);
        y_tmp = y .* f_i(x_tmp) ./ (1 - x_tmp .^ 2) .^ (0.5);
        c_i = 2 ./ pi .* trapz(x_tmp, y_tmp);
        y_i = y_i + c_i .* f_i(x_tmp);
        plot(x, y, x, y_i);
        legend('f(x)', 'i-ая частичная сумма ряда');
        xlabel('x');
        ylabel('y');
        title('Система функций Чебышёва');
        mov(i) = getframe();
    end
end
if strcmp(meth, 'Legendre')
    xBase = linspace(-1, 1, nDots);
    y_i = 0.5 * trapz(xBase, y);
    for i = 1:n
        f_i = getFuncLeg(i);
        y_tmp = y .* f_i(xBase);
        c_i = (2 * i + 1) ./ 2 .* trapz(xBase, y_tmp); 
        y_i = y_i + c_i .* f_i(xBase);
        plot(x, y, x, y_i);
        legend('f(x)', 'i-ая частичная сумма ряда');
        xlabel('x');
        ylabel('y');
        title('Система функций Лежандра');
        mov(i) = getframe();
    end
end
end