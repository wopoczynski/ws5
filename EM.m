function [out] = EM(sample)

%step 1 wszystkie mozliwe kobinacje SNP
[unphased,~,ic] = unique(sample, 'rows');
patientsAmount = accumarray(ic, 1);
result = [unphased patientsAmount];

%step 2 unikalne genotypy
tmp=[];
for i = 1 : size(result, 1) 
tmp{i} = translate2(result(i,1:3));
end

%step 3 unikalne haplotypy
genotype = char(tmp);
genotypeUnique = unique(genotype, 'rows');
index = 1:1:8;%indeksy haplotypow, potrzebuje ich numery
haplotypeStep3 = table(index', genotypeUnique);%tabelka slajd 38

%step 4 losowanie liczb sumuj¹cych siê do 1
randomNumber = rand(8,1);
for i=1:8%warunki pocz¹tkowe
initialFrequency(i) = randomNumber(i) / sum(randomNumber);
end
%tabelka slajd 39
haplotypeStep4 = table(index', genotypeUnique, initialFrequency)



%step 5 dzieje sie magia jeszcze nie wiem o co chodzi ale bêdzie pêtla w
%pêtli 


end
