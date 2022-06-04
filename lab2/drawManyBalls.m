function drawManyBalls(alphas, colors, edges)
lb = -3;
ub = 3;
nDots = 21;
x = linspace(lb, ub, nDots);
y = linspace(lb, ub, nDots);
z = linspace(lb, ub, nDots);
[X, Y, Z] = meshgrid(x, y, z);
if ((length(alphas) ~= length(colors)) || (length(alphas) ~= length(edges)))
    error('Размеры не совпадают');
end

for i = 1:length(alphas)
    figure;
    if isinf(alphas(i))
        V = max(max(abs(X), abs(Y)), abs(Z));
    else
        V = abs(X).^(alphas(i)) + abs(Y).^(alphas(i)) + abs(Z).^(alphas(i));
    end
    
    ball = isosurface(X,Y,Z,V,2);
    
    if (size(ball.vertices) == 0)
        error('Таких точек нет!');
    end
    paintedBall = patch(ball);
    xlabel('x');
    ylabel('y');
	zlabel('z');
    title('\alpha = ' + string(alphas(i)));
    grid on;
    paintedBall.FaceColor = colors(i);
    if (edges(i) == "None")
        paintedBall.EdgeColor = "None";
    else
        paintedBall.EdgeColor = edges(i);
    end
    daspect([1, 1, 1])
    view(3); 
    axis tight
    camlight('left');
end
end