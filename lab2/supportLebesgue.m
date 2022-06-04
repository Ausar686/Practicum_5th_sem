function [func] = supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub)
%f(x) - �������, ���������� ������
%A, Aeq - ������� 
%b, beq - ������� 
%c, ceq - �������, ������������ ������
%x0, lb, ub - ������� ��� �������
%������ min(f(x)): c(x) <= 0, ceq(x) = 0, A*x <= b, Aeq * x = beq, lb <= x <= ub
%���������� ������� ���� ������� �� t (������-�����������):
%2-�� �������� = max(<x, t>) = min(-<x, t>)
%1-�� �������� - �����, � ������� ����������� min(-<x, t>)
scalarProduct = @(x, t) x(1) * t(1) + x(2) * t(2);
function [c, ceq] = nlcon(x)
%����� ���������� ������� c � ceq (�� ���������� ����)
    c = f(x); %��������� ������� f(x) <= 0, ������� c(x) = f(x)
    ceq = [];
end
nonlcon = @nlcon;
fmin = @(t) fmincon(@(x) - scalarProduct(x,t), x0, A, b, Aeq, beq, lb, ub, nonlcon);
func = @(t) [fmin(t), scalarProduct(fmin(t),t)];
end