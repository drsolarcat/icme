%%%
% Noise-robust second derivative
%
% Alexey Isavnin
%%%

function ddy = NoiseRobustSecondDerivative(y, dx, N)
    
    
    M = (N-1)/2;
    s = zeros(1, M+2);
    if size(y, 2) < size(y, 1)
        s = s';
    end
    
    for i = M+1:-1:1
        k = i-1;
        if k == M
            s(i) = 1;
        else
            s(i) = ((2*N-10)*s(i+1)-(N+2*k+3)*s(i+2))/(N-2*k-1);
        end
    end

    ddy = zeros(size(y));
    ddy_f = zeros(size(y));
    ddy_b = zeros(size(y));
    for i = 1:length(ddy)
        if i == 1
            ddy(1) = (2*y(1)-5*y(2)+4*y(3)-y(4))/dx^2;
        elseif i == length(ddy)
            ddy(i) = (2*y(end)-5*y(end-1)+4*y(end-2)-y(end-3))/dx^2;
        elseif i > 1 && i <= M || i >= length(ddy)-M+1 && i < length(ddy)
            ddy(i) = (y(i+1)-2*y(i)+y(i-1))/dx^2;
        else
            ddy(i) = 1/2^(N-3)/dx^2*(s(1)*y(i)+sum(s(2:M+1).*(y(i+1:i+M)+y(i-1:-1:i-M))));
        end
    end 
    
    %{
    for i = 1:length(ddy)
        % forward
        if i <= length(ddy)-M-M
            ddy_f(i) = 1/2^(N-3)/dx^2*(s(1)*y(i+M)+sum(s(2:M+1).*(y(i+M+1:i+M+M)+y(i+M-1:-1:i))));
        end
        % backward
        if i >= M+M+1
            ddy_b(i) = 1/2^(N-3)/dx^2*(s(1)*y(i-M)+sum(s(2:M+1).*(y(i-M+1:i)+y(i-M-1:-1:i-M-M))));
        end
    end
    
    ddy_f(end-M-M+1:end) = ddy_b(end-M-M+1:end);
    ddy_b(1:M+M) = ddy_f(1:M+M);
    ddy = mean([ddy_f; ddy_b]);
    %}
end
