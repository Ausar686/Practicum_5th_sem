function res = func_1_7(a, b)
% ���������� ����������� ������ ��������:  abs(a_i - b_j)
res = max([max(a) - min(b), max(b) - min(a)]); % ������� ������� ����� ������� �������� ��������
end