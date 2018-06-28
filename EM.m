function out = EM(sample)

%step 1 wszystkie mozliwe kobinacje SNP
[unphased,~,ic] = unique(sample, 'rows');
patientsAmount = accumarray(ic, 1);
result = [unphased patientsAmount];

%step 2, 3 i 4 unikalne genotypy i haplotypy dla ka¿dego ze SNP
[possibleGenotype, code, initialFrequency] = translate2(result(:,1:3));
tableGenotype = possibleGenotypeTable(possibleGenotype, code);


%step 5 macierzP
matrixP = matrixPStep5(initialFrequency);

% step 6 prawdopodobieñstwa dla kazdego fenotypu
P = vectorP(tableGenotype, matrixP, unphased);

%step 7 macierz Pnew
matrixPnew = macierzPnewStep7(P, matrixP, tableGenotype);

%step 8
for i = 1 : size(code, 1)
    newEstimate(i,1) = 0.5 * (sum(matrixPnew(i,:), 2)+ matrixPnew(i,i));
end

%step 9 - wylicznie L i warunek stop

Lold = prod(P(:,1) .* P(:,2));
L=0; 
while abs(Lold - L)>0.01
matrixP = matrixPStep5(newEstimate);
P = vectorP(tableGenotype, matrixP, unphased);
matrixPnew = macierzPnewStep7(P, matrixP, tableGenotype);

for i = 1 : size(code, 1)
    newEstimate(i,1) = 0.5 * (sum(matrixPnew(i,:), 2)+ matrixPnew(i,i));
end
Lold = L;
L = prod(P(:,1) .* P(:,2));

end


nukleotydy = numToNuc(code);
out = table(nukleotydy, initialFrequency, newEstimate);


end

