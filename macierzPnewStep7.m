function matrixPnew = macierzPnewStep7(P, matrixP, tableGenotype)

N = sum(P(:,2));
for i = 1:size(matrixP,1)
    for j = 1:size(matrixP,1)
        idx = tableGenotype(:, 2) == i & tableGenotype(:, 3) == j;      
        if sum(idx) == 1
            numer = tableGenotype(idx,1);
            matrixPnew(i,j) = (P(numer,2)*matrixP(i,j))/(N*P(numer,1));
            matrixPnew(j,i) = (P(numer,2)*matrixP(i,j))/(N*P(numer,1));
        end
    end
end

end
