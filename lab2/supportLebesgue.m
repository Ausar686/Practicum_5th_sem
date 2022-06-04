function [func] = supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub)
%f(x) - функция, возвращает скаляр
%A, Aeq - матрицы 
%b, beq - векторы 
%c, ceq - функции, возвращающие вектор
%x0, lb, ub - векторы или матрицы
%Ищется min(f(x)): c(x) <= 0, ceq(x) = 0, A*x <= b, Aeq * x = beq, lb <= x <= ub
%Возвращает функцию пару функций от t (вектор-направление):
%2-ое значение = max(<x, t>) = min(-<x, t>)
%1-ое значение - точка, в которой достигается min(-<x, t>)
scalarProduct = @(x, t) x(1) * t(1) + x(2) * t(2);
function [c, ceq] = nlcon(x)
%Задаём нелинейные условия c и ceq (см примечание выше)
    c = f(x); %Выполняем условие f(x) <= 0, полагая c(x) = f(x)
    ceq = [];
end
nonlcon = @nlcon;
fmin = @(t) fmincon(@(x) - scalarProduct(x,t), x0, A, b, Aeq, beq, lb, ub, nonlcon);
func = @(t) [fmin(t), scalarProduct(fmin(t),t)];
end