function res = takeColumn(A,i)
    nStr = size(A,1);
    res = zeros(nStr,1);
    for j = 1:nStr
        res(j) = A(j,i);
    end
end