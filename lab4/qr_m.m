function [Q,R] = qr_m(A)
    column = @takeColumn; % берем столбец
    nCol = size(A,2); % количество столбцов
    b = zeros(nCol,nCol); % матрица столбцов
    c = zeros(1,nCol-1); % коэффициенты
    i = 1;
    rg = rank(A);
    if(norm(column(A,1)) > 10^(-5))
        b(:,1) = column(A,1)/norm(column(A,1));
    end
    Q = [b(:,1)];
    while(i < nCol)
            j = i+1;
            b(:,j) = column(A,j)';
            for k=1:i
                if(dot(column(b,k),column(b,k)) < 10^(-5)) % скалярное произведение
                    c(k) = 0;
                else
                    c(k) = dot(column(b,k),column(A,j))/dot(column(b,k),column(b,k));
                end
                b(:,j) = b(:,j) - c(k).*column(b,k);
            end
            if(norm(b(:,j)) > 10^(-5))
                b(:,j) = b(:,j)/norm(b(:,j));
            end
            i = i+1;
            Q = [Q b(:,j)];
    end
    
    if(Q == 0)
        R = -eye(nCol);
    else
        R = Q'*A;
    end
    
    if(rg < nCol)
        Q = eye(nCol);
        R = Q'*A;
    else
        Q = -Q;
        R = -R;
    end
end
