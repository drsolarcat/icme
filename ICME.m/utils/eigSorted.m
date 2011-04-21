%%%
% Get eigen values and vectors in some sorted order
%
% Alexey Isavnin
%%%

function [x1, x2, x3, lambda1, lambda2, lambda3] = eigSorted(M, order)
    [V, D] = eig(M);
    D = abs(D);
    D = [D(1, 1), D(2, 2), D(3, 3)];
    [D, order] = sort(D, order);
    % eigen vectors
    x1 = V(:, order(1));
    x2 = V(:, order(2));
    x3 = V(:, order(3));
    % eigen values
    lambda1 = D(1);
    lambda2 = D(2);
    lambda3 = D(3);
end
