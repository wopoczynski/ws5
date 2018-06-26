function XcSquare = HW(genAA, genAa, genaa, rows, chromosomAmount)
    p = (2 * genAA + genAa) / chromosomAmount; 
    q = 1 - p;
    pSquare = p.^2;
    pq2 = 2.*p.*q;
    qSquare = q.^2;
    observed = [genAA; genAa; genaa];
    observedHWE = [rows.*pSquare; rows.*pq2; rows.*qSquare];
    
    for i = 1:size(genAA,2)
        XcSquare(i) = (observed(1,i)-observedHWE(1,i))^2/observedHWE(1,i) + ...
            (observed(2,i)-observedHWE(2,i))^2/observedHWE(2,i) + ...
            (observed(3,i)-observedHWE(3,i))^2/observedHWE(3,i);
    end
