%dane 6
clear all;clc;close all;
%% Równowaga sprzê¿eñ 
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

load('Dane_Clark_EM/dane_6.mat');

[a] = clarkAlg(sample);

% algorytm em