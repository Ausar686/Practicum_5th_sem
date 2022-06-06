function res = func_1_7(a, b)
% Нахождение наибольшего модуля разности:  abs(a_i - b_j)
res = max([max(a) - min(b), max(b) - min(a)]); % Перебор случаев через крайние значения векторов
end
