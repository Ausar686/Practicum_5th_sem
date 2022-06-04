clc;
close all
GenerateTable(6, 8);
SaveTurtle (0);
clear;


function  SaveTurtle (var)
    minimum = @(action, vals) min(action + vals);
    % импортируем данные
    load('C:\Users\gleby\OneDrive\Документы\MATLAB\sem5\lab3\field','radiation', 'n', 'm');
    figure(1);
    ylim([0 m]); xlim([0 n]);
    xticks(1:n); yticks(1:m)
    hold on;
    
    % Рисуем поле
    Nmap = 10;
    cMin1 = [0 1 1];
    cMax1 = [1 1 0];
    cMap1 = zeros(Nmap, 3);
    for i=1:Nmap
        cMap1(i,:) = cMin1*(Nmap - i)/(Nmap - 1) + cMax1*(i - 1)/(Nmap - 1);
    end 
    for x = 1:n
        for y = 1:m
            fill([x-1 x-1 x x],[y-1 y y y-1], [radiation(x,y)/10, 1, 1-radiation(x,y)/10]);
        end
    end
    colormap(cMap1);
    colorbar('Ticks',0.05:0.1:1, 'TickLabels',{'1','2','3','4','5', '6','7','8','9','10'})
    
    % таблички с весами
    radiation = radiation+1;
    V = 1*ones(n+2,m+2)*inf;
    up = 1*ones(n+2,m+2); 
    down = 1*ones(n+2,m+2);
    left = 1*ones(n+2,m+2);
    right = 1*ones(n+2,m+2);
    % диагонали
    if var
        up_left = 1*ones(n+2,m+2); 
        up_right = 1*ones(n+2,m+2);
        down_left = 1*ones(n+2,m+2);
        down_right = 1*ones(n+2,m+2); 
    end
    
    for i = 1:n % left-right
        for j = 1:m %up-down
            up(i+1, j) = radiation(i, j);
            down(i+1, j+2) = radiation(i, j);
            left(i+2, j+1) = radiation(i, j);
            right(i, j+1) = radiation(i, j);
            
            % диагонали
            if var
                up_right(i,j) = radiation(i, j);
                up_left(i+2,j) = radiation(i, j);
                down_right(i,j+2) = radiation(i, j);
                down_left(i+2,j+2) = radiation(i, j);
            end 
        end
    end
    
    left([1 2 n+2],:) = inf;
    right([1 n+1 n+2],:) = inf;
    down(:,[1 2 m+2]) = inf;
    up(:,[1 m+1 m+2]) = inf;
    
    left = left.*1.5;
    down = down.*1.5;
    
    if var
        up_right = up_right.*0.8;
        down_left = down_left.*2;
        
        up_left([1 2 n+2],:) = inf;
        up_left(:,[1 m+1 m+2]) = inf;
        
        up_right([1 n+1 n+2],:) = inf;
        up_right(:,[1 m+1 m+2]) = inf;
        
        down_left(:,[1 2 m+2]) = inf;
        down_left([1 2 n+2],:) = inf;
        
        down_right(:,[1 2 m+2]) = inf;
        down_right([1 n+1 n+2],:) = inf;
    end
    
    V(n+1,m+1) = 0;
    n=n+1;
    m = m+1;
    iter_num = 1;
	change = 1*ones(n+1,m+1); % cells to change
    
    while (size(find(change == 1), 1) > m*2 + n*2 + 1)
        for i = n:-1:2
            for j = m:-1:2
                if ((j ~= m) || (i ~= n))&&(change(i, j) == 1)
                    if var
                        val = minimum([right(i,j) left(i, j) up(i, j) down(i, j) ...
                                       up_right(i,j) up_left(i, j) down_right(i, j) down_left(i, j)], ...
                                      [V(i+1, j)  V(i-1, j) V(i, j+1) V(i, j-1) ...
                                       V(i+1, j+1)  V(i-1, j+1) V(i+1, j-1) V(i-1, j-1)]);
                        if V(i, j) ~= val
                            change([i-1 i-1 i+1 i+1] , [j-1 j+1 j-1 j+1]) = 1;
                            change([i-1 i+1 i i] , [j j j-1 j+1]) = 1;
                        end
                        V(i, j) = val;
                        
                    else    
                        val = minimum([right(i,j) left(i, j) up(i, j) down(i, j)], ...
                                      [V(i+1, j)  V(i-1, j) V(i, j+1) V(i, j-1)]);
                        if V(i, j) ~= val
                            change([i-1 i+1 i i] , [j j j-1 j+1]) = 1;
                        end
                        V(i, j) = val;
                    end 
                    change(i, j) = 0;
                end
            end
        end
        iter_num = 1 + iter_num;
        size(find(change == 1))
    end
    str = strcat('Количество итераций = ', num2str(iter_num));
    disp(str)
    
    i = 2;
    j = 2;
    k = 0;
    total_val = 0;
    while (k < n*m)
        k = k+1;
        if (i == n)&&(j == m)
            str = strcat('Количество шагов = ', num2str(k-1));
            title (str);
            disp(str)
            break;
        end
        
        if var
            [val, ind] = minimum([right(i,j) left(i, j) up(i, j) down(i, j) ...
                                 up_right(i,j) up_left(i, j) down_right(i, j) down_left(i, j)], ...
                                 [V(i+1, j)  V(i-1, j) V(i, j+1) V(i, j-1) ...
                                  V(i+1, j+1)  V(i-1, j+1) V(i+1, j-1) V(i-1, j-1)]);
        else
            [val, ind] = minimum([right(i,j) left(i, j) up(i, j) down(i, j)], ...
                                [V(i+1, j)  V(i-1, j) V(i, j+1) V(i, j-1)]);
        end
        total_val = total_val + val;
        plot([i i+1*(ind==1)-1*(ind==2)+1*(ind==5)-1*(ind==6)+1*(ind==7)-1*(ind==8)]-1.5, ...
             [j j+1*(ind==3)-1*(ind==4)+1*(ind==5)+1*(ind==6)-1*(ind==7)-1*(ind==8)]-1.5, 'b*-' )
        i = i+1*(ind==1)-1*(ind==2)+1*(ind==5)-1*(ind==6)+1*(ind==7)-1*(ind==8);
        j = j+1*(ind==3)-1*(ind==4)+1*(ind==5)+1*(ind==6)-1*(ind==7)-1*(ind==8);
        
    end
    hold off
end

function radiation = GenerateTable (n, m)
    radiation  = round(rand(n, m).*10);
    %radiation  = round(ones(n, m).*0); %везде одинаковая радиация
    save('C:\Users\gleby\OneDrive\Документы\MATLAB\sem5\lab3\field','radiation', 'n', 'm');
end
