function matrixP = matrixPStep5(initialFrequency)
matrixP = zeros(size(initialFrequency,1), size(initialFrequency,1));
for i=1:size(initialFrequency,1)
    for j=1:size(initialFrequency,1)
        if (j == i)
            matrixP(i,j) = initialFrequency(i)^2;
        else
            matrixP(i,j) = 2*initialFrequency(i)*initialFrequency(j);
        end
    end
end
end
