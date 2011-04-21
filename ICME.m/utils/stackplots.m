%%%
%
% Place subplots in a pretty vertical stack with shared x axis
%
% Alexey Isavnin
%
%%%

function stackplots(fhandle, d)

    fpos = get(fhandle, 'Position');

    ahandles = get(fhandle, 'children');
    N = length(ahandles);
    set(ahandles, 'Units', 'normalized');
    set(ahandles, 'TickDir', 'out');
    set(ahandles, 'Box', 'on');
    
    d = d/fpos(4);
    
    xaxis = copyobj(ahandles(1), fhandle);
    delete(get(xaxis, 'children'));
    delete(get(xaxis, 'YLabel'));
    set(xaxis, 'YColor', get(fhandle, 'Color'));
    set(xaxis, 'Color', get(fhandle, 'Color'));
    set(xaxis, 'Box', 'off');
    
    
    set(ahandles(1:N), 'XTick', []);
    for i = 1:N
        set(ahandles(i), 'XTickLabel', []);
        set(ahandles(i), 'XColor', get(ahandles(i), 'Color'));
        delete(get(ahandles(i), 'XLabel'));
        delete(get(ahandles(i), 'Title'));
    end
    
    pxt = get(xaxis, 'TightInset');
    pxo = get(xaxis, 'OuterPosition');
    pxo = [pxo(1) 0 pxo(3) d+pxt(2)+pxt(4)];
    set(xaxis, 'OuterPosition', pxo);
    h = (1-pxo(2)-pxo(4)-d*N)/N;
    xlim = get(xaxis, 'XLim');
    ylblx = xlim(1);
    for i = 1:N
        p = get(ahandles(i), 'Position');
        set(ahandles(i), 'Position', [p(1) pxo(4)+h*(i-1)+d*(i-1) p(3) h]);
        ylbl = get(ahandles(i), 'YLabel');
        p = get(ylbl, 'Position');
        if p(1) < ylblx
            ylblx = p(1);
        end
    end
    
    for i = 1:N
        ylbl = get(ahandles(i), 'YLabel');
        p = get(ylbl, 'Position');
        p(1) = ylblx;
        set(ylbl, 'Position', p);
    end
    
end
