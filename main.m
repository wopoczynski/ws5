%dane 6
clear all;clc;close all;
%% R�wnowaga sprz�e� 
load('Dane_LD/dane_6_LD.mat');

[counts, frequency, pq, recombination, P] = allelFrequency(sample);

[chi2, D, chromosomCorelation] = X2(pq, P);
%XcSquare [snp1;snp2;snp3];

% chyba tak to trzeba policzyc dla kazdego z snp
XSquareHealthy = HW(counts(:,:,1));
XSquareSic = HW(counts(:,:,2));
%XcSquare [snp1;snp2;snp3];

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
% [haplotypeHealthy, freqHealthy] = clarkAlg(healthy);
% [haplotypeSic, freqSic] = clarkAlg(healthy);
% 
% [translatedH] = numToNuc(haplotypeHealthy)
% [translatedS] = numToNuc(haplotypeSic)

% algorytm em

healthyEM = EM(healthy)
sicEM = EM(sic)

%todo

a = (141.25+76.25+13.5+2.75+2.25+27+54.75+57.25)


