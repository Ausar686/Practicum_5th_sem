function  drawSet(rho, N)
    figure;
    alpha = linspace(0, 2 * pi, N);
    xIn = zeros(1, N);            
    yIn = zeros(1, N);  
    xOut = zeros(1, N);
    yOut = zeros(1, N);
    %Инициализация векторов для вычисление опорной функции
    v1 = [0, 0]; 
    v2 = [0, 0]; 
    for i = 1:N
        v1 = v2;              
        v2 = [cos(alpha(i)), sin(alpha(i))]; 
        z = rho(v2);      
        xIn(i) = z(1);
        yIn(i) = z(2);

        if (i == 1)
            c2 = -v2(1)*xIn(i) -v2(2)*yIn(i); 
        else
            c1 = c2;                
            c2 = -v2(1)*xIn(i) - v2(2)*yIn(i);
            %Магические формулы для точки пересечения касательных
            yOut(i - 1) =(v1(1) * c2 - v2(1) * c1) / (-v1(1) * v2(2) + v1(2) * v2(1));
            xOut(i - 1) =(v1(2) * c2 - v2(2) * c1) / (-v1(1) * v2(2) + v1(2) * v2(1));
        end
    end
    yOut(N) = yOut(1);
    xOut(N) = xOut(1);
    plot(xIn, yIn, '--', xOut, yOut, '.-');
    axis tight;
    grid on;
    title(func2str(rho));
    xlabel('x');
    ylabel('y');
    legend('Внутренняя аппроксимация', 'Внешняя аппроксимация');
end