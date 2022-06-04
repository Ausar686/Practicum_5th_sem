%close all;
clc;
func1 = @(t) 2.*t./(t.^2 - 4.*t + 5);

func3 = @(t) t.^3.*exp(-t.^2 - t);
func4 = @(t) cos(t)./(1 + 4.*abs(t).^5);


%close all
fighandle1 = figure('Name','Function 1, empty axes','NumberTitle','off');
res = plotFT(fighandle1, func1, @(t) ftfunc1(t) , 0.01, [-10 10], [-30 30]);
%pause

%close all
fighandle1 = figure('Name','Function 2, empty axes','NumberTitle','off');
res = plotFT(fighandle1,  @(t) func2(t),  @(t) ftfunc2(t) , 0.01, [-10 10], [-50 50]);
%pause


%close all
fighandle1 = figure('Name','Function 3, empty axes','NumberTitle','off');
res = plotFT(fighandle1, func3, [], 1, [-10 10], [-50 50]);
%pause

%close all
fighandle1 = figure('Name','Function 4, empty axes','NumberTitle','off');
res = plotFT(fighandle1, func4, [] , 1, [-10 10], [-50 50]);
%pause

clear;





function res = plotFT(hFigure, fHandle, fFTHandle, step, inpLimVec, outLimVec)
    figure(hFigure); % make a current figure
    if isfield(hFigure.UserData, 'inpLimVec')
        a = hFigure.UserData.inpLimVec(1);
        b = hFigure.UserData.inpLimVec(2);
    else
        a = inpLimVec(1);
        b = inpLimVec(2);
    end
    
    if isfield(hFigure.UserData, 'outLimVec')
        c = hFigure.UserData.outLimVec(1);
        d = hFigure.UserData.outLimVec(2);
    elseif length(outLimVec) > 1 
        c = outLimVec(1);
        d = outLimVec(2);
    elseif length(hFigure.Children) > 1
        c = hFigure.Children.XLim(1);
        d = hFigure.Children.XLim(2);
    else
        % частота Найквиста
        c = -1/(2*step);
        d = 1/(2*step); 
    end
    
    if length(hFigure.Children) > 1
        real_ax = hFigure.Children(1);
        im_ax = hFigure.Children(2);
    else
        real_ax = axes('Parent', hFigure, 'OuterPosition',[0    0.5307    1.0000    0.4248]);
        im_ax = axes('Parent', hFigure, 'OuterPosition',[0    0.0568    1.0000    0.4146]);
    end
    real_ax.NextPlot = 'add';
    im_ax.NextPlot = 'add';
    SPlotInfo = get(hFigure, 'UserData');
    if ~isempty(SPlotInfo) && isfield(hFigure.UserData, 'ax_handles')
        delete(SPlotInfo.ax_handles);
    end
    
    T = b - a;  
    N = T/step;
    N = round(N);
    %N = 2^nextpow2(N)
    step = T / N;
    x = a:step:b;
    y = fHandle(x);
    if ~isempty(find(isnan(y), 1))
        find(isnan(y))
        y(isnan(y)) = 0;
    end
    n = 0;
    if (sign(a) + sign(b) == 2)
        n = round(a/T);
    elseif (sign(a) + sign(b) == -2)
        n = -round(b/T);
    end
    border = n*T;
    ind = find(x <= border, 1, 'last' );
    Y_new(1:N+1-ind) = y(ind+1:N+1);
    Y_new(N+1-ind+1:N+1) = y(1:ind);
    Yft = step.*fft(Y_new);%.*exp(1i.*(b/2 - a/2).*linspace(a, b, size(Y_new, 1)));
    
    stepft = 2*pi/T;
    Xft = 0:stepft:stepft*N;

    hold on
    Tft = stepft*N; 
    lborder = -Tft;
    rborder = Tft; 
    counter = 2;
    while (c < lborder || d > rborder)
        lborder = lborder - Tft;
        rborder = rborder + Tft;
        counter = counter + 2; 
    end
    
    
    Xft_new = lborder+stepft:stepft:rborder;
    Yft_new = repmat(Yft(2:end), 1, counter);
    re_y = real(Yft_new);
    im_y = imag(Yft_new);
    
    plot(Xft_new, re_y, 'r',  'Parent', real_ax);
    plot(Xft_new, im_y, 'r', 'Parent', im_ax);

    if isempty(fFTHandle) == 0
        xf = c:0.1:d;
        re_f = real(fFTHandle(xf));
        im_f = imag(fFTHandle(xf));
        plot (xf, re_f, 'b', 'Parent', real_ax)
        plot (xf, im_f, 'b', 'Parent', im_ax)
        real_ax.YLim = [-max([re_y re_f]) max([re_y re_f])*1.2];
        im_ax.YLim = [min([im_y im_f])*1.2 max([im_y im_f])*1.2];
        legend ('Численное решение','Аналитическое решение');
    else 
        real_ax.YLim = [min(re_y) - abs(min(re_y)) max(re_y)];
        im_ax.YLim = [min(im_y) - abs(min(im_y)) max(im_y)];
        legend ('Численное решение')
    end
    real_ax.XLim = [c d];
    im_ax.XLim = [c d];
    
    real_ax.XLabel.String = '\lambda';
    im_ax.XLabel.String = '\lambda';
    im_ax.YLabel.String = 'Im(F(\lambda))';
    real_ax.YLabel.String = 'Re(F(\lambda))';
    real_ax.Title.String = 'Вещественная часть';
    im_ax.Title.String =  'Мнимая часть';
    SPlotInfo = struct('ax_handles', [real_ax im_ax], 'outLimVec', [a b], 'inpLimVec', [c d], 'nPoints', n);
    hold off;
    set(hFigure,'UserData',SPlotInfo);
    res.inpLimVec = [a b];
    res.outLimVec = [c d];
end

function f = ftfunc1 (t)
    f = (2 .* pi) .* exp((-1 - 2.*i).*t) .* (((2 + i) .* exp(2.*t).*(t<0)) + (2 - i) .*(t>0));
    % ПФ ф-ии func1
end

function f = ftfunc3 (t)
    f = 1./64.*sqrt(pi/2).*exp(-(t - i).^2./8).*(1 + i.*t).*(-13+t.*(t - 2.*i));
    % ПФ ф-ии func1
end


function f = func2 (t)
    f = t.^2 + t;
    f(abs(t) >= 1) = 0;
end

function f = ftfunc2 (t)
    f = (exp(i.*t) - exp(-i.*t))./(i.*t) + (1 + 2./(i.*t)).*(exp(-i.*t).*(i.*t +1) - exp(i.*t).*(1-i.*t))./t.^2;
    % ПФ ф-ии func2
end

function f = dirac (t)
    f = t*0;
    f(abs(t) == 0) = inf;
end
