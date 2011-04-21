%%%
%
% vline across multiple subplots
%
% Alexey Isavnin
%
%%%

function h = subplotvline(x, in1, in2)

    if nargin < 3
        in2 = [];
        if nargin == 1
            in1 = [];
        end
    end
    
    ahandles = get(gcf, 'children');
    xlim = get(ahandles(end), 'XLim');
    xx = (x-xlim(1))./(xlim(2)-xlim(1));
    units = get(gcf, 'units');
    set(gcf, 'units', 'normalized');
    miny = 1;
    maxy = 0;
    minx = 1;
    maxx = 0;
    for k = 1:length(ahandles)
        p = get(ahandles(k), 'Position');
        if p(2) < miny
            miny = p(2);
        end
        if (p(2)+p(4)) > maxy
            maxy = p(2)+p(4);
        end
        if p(1) < minx
            minx = p(1);
        end
        if (p(1)+p(3)) > maxx
            maxx = p(1)+p(3);
        end
    end
    
    axes('visible', 'off', 'position', [minx miny maxx-minx maxy-miny]);
    axis manual;
    h = vline(xx, in1, in2);
    set(gcf, 'units', units);

end
