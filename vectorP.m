function P = vectorP(tableGenotype, matrixP, unphased)


P = zeros(size(unphased,1),1);
for i = 1:size(unphased,1)
    idx = tableGenotype(:,1) == i;
    x = tableGenotype(idx,2:3);
    for j = 1:size(x,1)
        idx1 = x(j,1);
        idx2 = x(j,2);
        P(i, 1) = P(i, 1) + matrixP(idx1, idx2);
        P(i,2) = sum(idx);
    end
end

end