function res = func_1_7(a, b)
% Ќахождение наибольшего модул€ разности:  abs(a_i - b_j)
res = max([max(a) - min(b), max(b) - min(a)]); % ѕеребор случаев через крайние значени€ векторов
end