function plotRM()
    theta = 0:1:90;
    phi = 0:1:360;
    plotRM_('../Science/icme/ICME.cpp/rm_original.dat', theta, phi);
    plotRM_('../Science/icme/ICME.cpp/rm_combined.dat', theta, phi);
end

function plotRM_(path, theta, phi)
    r = load(path);
    R = r.^-1;
    figure
    h = polar([0 2*pi], [0 90]);
    ph = findall(gca, 'type', 'patch');
    set(ph, 'facecolor', [0 0 .5], 'edgecolor', [0 0 .5]);
    set(gcf, 'Color', 'white');
    set(gcf, 'Renderer', 'Painters');
    pl = findobj(allchild(gca));
    hold on
    contour(theta'*cos(phi*pi/180), theta'*sin(phi*pi/180), R, ...
            linspace(min(R(:)), max(R(:)), 500), 'Fill', 'on');
    colorbar
    for i = 1:length(pl)-1
        if strcmpi(get(pl(i), 'Type'), 'line')
            set(pl(i), 'Color', 'white');
        elseif strcmpi(get(pl(i), 'Type'), 'text') && i > 25
            set(pl(i), 'Color', 'white');
        end

        uistack(pl(i), 'top');
    end
    delete(h)
    cbar_handle = findobj(gcf, 'Tag', 'Colorbar');
    set(get(cbar_handle,'xlabel'),'string','1/R');
    hold off
end