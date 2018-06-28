%dane 6
clear all;clc;close all;
%% R�wnowaga sprz�e� 
load('Dane_LD/dane_6_LD.mat');

[counts, frequency, pq, recombination, P] = allelFrequency(sample);

[chi2, D, chromosomCorelation] = X2(pq, P);

% chyba tak to trzeba policzyc dla kazdego z snp
XcSquare = zeros(3,1);
for i = 1:3
    %               AA              Aa              aa             rows amount              chromosomy
    XcSquare(i,:) = HW(counts(i,1,1), counts(i,2,1), counts(i,2,1), sum(counts(1,:,1)), sum(counts(1,:,1)) *2 );
end

chi2Critical = 5.9915;
%%
% algorytm clarka
clc;close all;clear all;
load('Dane_Clark_EM/dane_6.mat');

tmp = sample(sample(:,1) == 1,:);
healthy = tmp(:,3:end);
tmp = sample(sample(:,1) == 0,:);
sic = tmp(:,3:end);  

%% nie puszcza� ca�o�ci tylko pojedy�czo i najlepiej kilka razy - permutacje tworz�ce losowo�� czasem powoduj� wypadanie z tablic a to matlab i jest koniec 
[haplotypeHealthy, freqHealthy] = clarkAlg(healthy);
[haplotypeSic, freqSic] = clarkAlg(healthy);

[translatedH] = numToNuc(haplotypeHealthy)
[translatedS] = numToNuc(haplotypeSic)

% algorytm em

%todo



