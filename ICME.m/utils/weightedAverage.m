%%%
% Apply weighted 3point running average to some vector of data, where w is weight function
% values and N is number of runs to perform
%
% Alexey Isavnin
%%%

function A = weightedAverage(A, w, N)
    for i = 1:N
        A(2:end-1) = w*A(2:end-1)+1/2*(1-w)*(A(1:end-2)+A(3:end));
    end
end
