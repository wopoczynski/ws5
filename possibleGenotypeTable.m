function [table] = possibleGenotypeTable(possibleGenotype, code)

for i = 1:2:size(possibleGenotype,1)
    for j = 1:size(code,1)
        if possibleGenotype(i,1:3) == code(j,:)
            table(i,2) = j;%hap1
        end
        if possibleGenotype(i+1,1:3)== code(j,:)
            table(i,3) = j;%hap2
        end
        table(i,1) = (possibleGenotype(i,4));%fenotyp
    end
end
idx = sum(table,2) == 0;
table(idx,:)= [];

end