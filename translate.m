function [knownHaplo] = translate(knownHaplotypes, knownHaplo)

    knownHaploTmp = [];
    knownHaploTmp2 = [];
    
    snp10 = 'G';
    snp11 = 'A';
    snp11a = 'G';
    snp12 = 'A';

    snp20 = 'C';
    snp21 = 'C';
    snp21a = 'T';
    snp22 = 'T';

    snp30 = 'G';
    snp31 = 'A';
    snp31a = 'G';
    snp32 = 'A';
    
    %snp1
    if knownHaplotypes(1,1) == 0
        knownHaploTmp = [knownHaploTmp snp10];
    end
    if knownHaplotypes(1,1) == 1
        knownHaploTmp = [knownHaploTmp snp11];
        knownHaploTmp2 = [knownHaploTmp2 snp11a];
    end
    if knownHaplotypes(1,1) == 2
        knownHaploTmp = [knownHaploTmp snp12];
    end
    %snp2
    if knownHaplotypes(1,2) == 0
        knownHaploTmp = [knownHaploTmp snp20];
         if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp32];
        end
    end
    if knownHaplotypes(1,2) == 1
        knownHaploTmp = [knownHaploTmp snp21];
        knownHaploTmp2 = [knownHaploTmp2 snp21a];
    end
    if knownHaplotypes(1,2) == 2
        knownHaploTmp = [knownHaploTmp snp22];
         if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp32];
        end
    end
    
    %snp3
    if knownHaplotypes(1,3) == 0
        knownHaploTmp = [knownHaploTmp snp30];
         if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp32];
        end
    elseif knownHaplotypes(1,3) == 1
        knownHaploTmp = [knownHaploTmp snp31];
        knownHaploTmp2 = [knownHaploTmp2 snp31a];
    elseif knownHaplotypes(1,3) == 2
        knownHaploTmp = [knownHaploTmp snp32];
        if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp32];
        end
    end
    
    %sprawdzenie czy jakiœ nowy jest
    if (~any(strcmp(knownHaplo,knownHaploTmp)))
        knownHaplo{1,end+1} = knownHaploTmp;
        if numel(knownHaploTmp2)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
    end
end