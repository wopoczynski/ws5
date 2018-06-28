function [translated] = numToNuc(haplotypeHealthy)
translated = {};
G = 3;
A = 1;
C = 2;
T = 4;
for i = 1:size(haplotypeHealthy,1)
    for j = 1:size(haplotypeHealthy,2)
        if haplotypeHealthy(i,j) == A
            translated(i,j) = {'A'};
        elseif haplotypeHealthy(i,j) == C
            translated(i,j) = {'C'};
        elseif haplotypeHealthy(i,j) == G
            translated(i,j) = {'G'};
        elseif haplotypeHealthy(i,j) == T
            translated(i,j) = {'T'};
        end
    end
end
