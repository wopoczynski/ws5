

load('Dane_LD/dane_6_LD.mat');

%step 1 wszystkie mozliwe kobinacje SNP
[unphased,~,ic] = unique(sample(:,3:5), 'rows');
patientsAmount = accumarray(ic, 1);
result = [unphased patientsAmount];

%step 2, 3 i 4 unikalne genotypy i haplotypy dla ka¿dego ze SNP
[possibleGenotype, code, initialFrequency] = translate2(result(:,1:3));
tableGenotype = possibleGenotypeTable(possibleGenotype, code);


%step 5 macierzP
tableMatrixP = zeros(size(initialFrequency,1), size(initialFrequency,1));
for i=1:size(initialFrequency,1)
    for j=1:size(initialFrequency,1)
        if (j == i)
            tableMatrixP(i,j) = initialFrequency(i)^2;
        else
            tableMatrixP(i,j) = 2*initialFrequency(i)*initialFrequency(j);
        end
    end
end

% step 6 prawdopodobieñstwa dla kazdego fenotypu
Pi = zeros(size(unphased,1),1);
for i = 1:size(unphased,1)
    idx = tableGenotype(:,1) == i;
    x = tableGenotype(idx,2:3);
    for j = 1:size(x,1)
        idx1 = x(j,1);
        idx2 = x(j,2);
        Pi(i) = Pi(i) + tableMatrixP(idx1, idx2);
    end
end

%step 7

N = sum(patientsAmount);
for i = 1:size(code,1)
    for j = 1:size(code,1)
        idx = tableGenotype(:, 2) == i & tableGenotype(:, 3) == j;
        if sum(idx) == 1
            number = tableGenotype(idx,1);
            Pnew(i,j) = (patientsAmount(number)*tableMatrixP(i,j))...
                /(N * Pi(number));
        else
            Pnew(i,j) = 0;
        end
    end
end

%step 8
for i = 1 : size(code, 1)
    newEstimate(i,1) = 0.5 * (sum(Pnew(i,:), 2)+ Pnew(i,i));
end

%step 9 


