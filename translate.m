function [knownHaplo] = translate(knownHaplotypes, knownHaplo)

    knownHaploTmp = [];
    knownHaploTmp2 = [];
    knownHaploTmp3 = [];
    knownHaploTmp4 = [];
    knownHaploTmp5 = [];
    knownHaploTmp6 = [];
    knownHaploTmp7 = [];
    knownHaploTmp8 = [];

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
    elseif knownHaplotypes(1,1) == 1
        knownHaploTmp = [knownHaploTmp snp11];
        knownHaploTmp2 = [knownHaploTmp2 snp11a];
    elseif knownHaplotypes(1,1) == 2
        knownHaploTmp = [knownHaploTmp snp12];
    end
    
    %snp2
    if knownHaplotypes(1,2) == 0
        knownHaploTmp = [knownHaploTmp snp20];
         if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp20];
        end
    elseif knownHaplotypes(1,2) == 1
        tmp = knownHaploTmp;
        tmp2 = knownHaploTmp2;
        knownHaploTmp = [tmp snp21];
        knownHaploTmp2 = [tmp2 snp21a];
        knownHaploTmp3 = [tmp snp21a];
        knownHaploTmp4 = [tmp2 snp21];
    elseif knownHaplotypes(1,2) == 2
        knownHaploTmp = [knownHaploTmp snp22];
         if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp22];
        end
    end
    
    %snp3
    if knownHaplotypes(1,3) == 0
        knownHaploTmp = [knownHaploTmp snp30];
         if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp30];
         end
         if numel(knownHaploTmp3)>0
            knownHaploTmp3 = [knownHaploTmp3 snp30];
         end
        if numel(knownHaploTmp4)>0
            knownHaploTmp4 = [knownHaploTmp4 snp30];
        end
    elseif knownHaplotypes(1,3) == 1
        tmp = knownHaploTmp;
        tmp2 = knownHaploTmp2;
        tmp3 = knownHaploTmp3;
        tmp4 = knownHaploTmp4;
        knownHaploTmp = [tmp snp31];
        knownHaploTmp5 = [tmp snp31a];
        if numel(knownHaploTmp3)>0
            knownHaploTmp3 = [tmp3 snp31];
            knownHaploTmp6 = [tmp3 snp31a];
         end
        if numel(knownHaploTmp4)>0
            knownHaploTmp4 = [tmp4 snp31];
            knownHaploTmp7 = [tmp4 snp31a];
        end
        knownHaploTmp2 = [tmp2 snp31];
        knownHaploTmp8 = [tmp2 snp31a];
    elseif knownHaplotypes(1,3) == 2
        knownHaploTmp = [knownHaploTmp snp32];
        if numel(knownHaploTmp2)>0
            knownHaploTmp2 = [knownHaploTmp2 snp32];
        end
        if numel(knownHaploTmp3)>0
            knownHaploTmp3 = [knownHaploTmp3 snp32];
         end
        if numel(knownHaploTmp4)>0
            knownHaploTmp4 = [knownHaploTmp4 snp32];
        end
    end
    
    %sprawdzenie czy jakiœ nowy jest
    if (~any(strcmp(knownHaplo,knownHaploTmp)))
        knownHaplo{1,end+1} = knownHaploTmp;
        if numel(knownHaploTmp2)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
        if numel(knownHaploTmp3)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
        if numel(knownHaploTmp4)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
        if numel(knownHaploTmp5)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
        if numel(knownHaploTmp6)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
        if numel(knownHaploTmp7)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
        if numel(knownHaploTmp8)>0
            knownHaplo{1,end+1} = knownHaploTmp2;
        end
    end
