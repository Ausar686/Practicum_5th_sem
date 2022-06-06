function [res] = getEqual(f, g, t0, t1, N)
figure;
if (t1 <= t0)
    error('T1 <= T0');
end
h = (t1 - t0) / (N - 1);
t = linspace(t0, t1, N); % равномерное разбиение
%matrix - матрица точек p
matrix(1:N, 1) = f(t);
matrix(1:N, 2) = g(t);
%allDist - все возможные расстояния между p_i и p_j
allDist = pdist2(matrix, matrix); 
%dist - расстояния от p_i до p_{i+1}- сдвинутая диагональ matrix
dist = diag(allDist,1);
%midDist - среднее расстояние от p_i до p_{i+1}
midDist = sum(dist) / (N - 1);
disp('Uniform average dist');
disp (midDist);
%newCoords - новые координаты для поиска искомых значений t_1,...t_N
newCoords = t;
i = 1;
j = 5 * N;
while (max(dist) - min(dist)) > 0.01
    w_j = (t1 - t0) / j;
    if dist(i) < midDist
        if i == N - 1
            if (newCoords(i) - w_j) > newCoords(i - 1)
                newCoords(i) = newCoords(i) - w_j;
                matrix(1:N, 1) = f(newCoords);
                matrix(1:N, 2) = g(newCoords);
                allDist = pdist2(matrix, matrix); 
                dist = diag(allDist,1);
                midDist = sum(dist) / (N - 1);
            end
        else
            if ((newCoords(i + 1) + w_j) < newCoords(i + 2)) && ((newCoords(i + 1) + w_j) > newCoords(i))
                newCoords(i + 1) = newCoords(i + 1) + w_j;
                matrix(1:N, 1) = f(newCoords);
                matrix(1:N, 2) = g(newCoords);
                allDist = pdist2(matrix, matrix); 
                dist = diag(allDist,1);
                midDist = sum(dist) / (N - 1);
            end
        end
    elseif dist(i) > midDist
        if i == N - 1
            if (newCoords(i) + w_j) < newCoords(i+1)
                newCoords(i) = newCoords(i) + w_j;
                matrix(1:N, 1) = f(newCoords);
                matrix(1:N, 2) = g(newCoords);
                allDist = pdist2(matrix, matrix); 
                dist = diag(allDist,1);
                midDist = sum(dist) / (N - 1);
            end
        else
            if (newCoords(i + 1) - w_j) > newCoords(i)
                newCoords(i + 1) = newCoords(i + 1) - w_j;
                matrix(1:N, 1) = f(newCoords);
                matrix(1:N, 2) = g(newCoords);
                allDist = pdist2(matrix, matrix); 
                dist = diag(allDist,1);
                midDist = sum(dist) / (N - 1);
            end
        end
    end
    if i == N - 1
        i = 0;
        %j = j + 1;
    end

    i = i + 1;
    j = j + 1;
end
res(1, 1:N) = f(newCoords);
res(2, 1:N) = g(newCoords);
nDots = 1001;
t = linspace(t0, t1, nDots);
plot(f(newCoords), g(newCoords), '-o', f(t), g(t), f(t(1)), g(t(1)), 'c*', f(t(nDots)), g(t(nDots)), 'g*'); 
xlabel('f(t');
ylabel('g(t)');
disp('Required distance is');
disp(mean(dist));
axis equal;
end
