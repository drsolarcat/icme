%%%
% Project a set of vectors to some coordinate system
%
% Alexey Isavnin
%%%

function [Px, Py, Pz] = projection(Ax, Ay, Az, x, y, z)

    if size(x, 1) ~= 3
        x = x';
    end
    
    if size(y, 1) ~= 3
        y = y';
    end
    
    if size(z, 1) ~= 3
        z = z';
    end
    
    if size(Ax, 2) > size(Ax, 1)
        Ax = Ax';
    end
    
    if size(Ay, 2) > size(Ay, 1)
        Ay = Ay';
    end
    
    if size(Az, 2) > size(Az, 1)
        Az = Az';
    end

    M = length(Ax);
    
    if size(x, 2) == 1
        x = x*ones(1, M);
    end
    
    if size(y, 2) == 1
        y = y*ones(1, M);
    end
    
    if size(z, 2) == 1
        z = z*ones(1, M);
    end
    
    A = [Ax, Ay, Az]';
    Px = dot(A, x)';
    Py = dot(A, y)';
    Pz = dot(A, z)';
end
