%%%
% Construction of the magnetic variance matrix
% 
% Alexey Isavnin
%%%

function M = getVarianceMatrix(Bx, By, Bz)
    M = zeros(3, 3);
    M(1, 1) = mean(Bx.*Bx) - mean(Bx)*mean(Bx);
    M(1, 2) = mean(Bx.*By) - mean(Bx)*mean(By);
    M(1, 3) = mean(Bx.*Bz) - mean(Bx)*mean(Bz);
    M(2, 1) = mean(By.*Bx) - mean(By)*mean(Bx);
    M(2, 2) = mean(By.*By) - mean(By)*mean(By);
    M(2, 3) = mean(By.*Bz) - mean(By)*mean(Bz);
    M(3, 1) = mean(Bz.*Bx) - mean(Bz)*mean(Bx);
    M(3, 2) = mean(Bz.*By) - mean(Bz)*mean(By);
    M(3, 3) = mean(Bz.*Bz) - mean(Bz)*mean(Bz);
end
