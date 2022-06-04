function drawBall(alpha, level, params)
figure;
x = linspace(params.ax, params.bx, params.nDots);
y = linspace(params.ay, params.by, params.nDots);
z = linspace(params.az, params.bz, params.nDots);
[X, Y, Z] = meshgrid(x, y, z);
if isinf(alpha)
    V = max(max(abs(X), abs(Y)), abs(Z));
else
    V = abs(X).^(alpha) + abs(Y).^(alpha) + abs(Z).^(alpha);
end
ball = isosurface(X,Y,Z,V,level);
if (size(ball.vertices) == 0)
    error('Таких точек нет!');
end

paintedBall = patch(ball);
xlabel('x');
ylabel('y');
zlabel('z');
title('Линия уровня');
grid on;
paintedBall.FaceColor = params.faceColor;
paintedBall.EdgeColor = params.edgeColor;
daspect([1, 1, 1])
view(3); 
axis tight
camlight('left');
end