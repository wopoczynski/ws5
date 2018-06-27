function [out] = clarkAlg(sample)

  %% krok 1 zliczanie wystêpuj¹cych i losowanie kolejnosci
    [unphased,~,ic] = unique(sample, 'rows');
    patientsAmount = accumarray(ic, 1);
    result = [unphased patientsAmount];
    n = size(result,1);
    idx = randperm(n);
    result = result(idx,:);

    unphased = result(:,1:end-1); % mozliwosci
    patientsAmount = result(:,end); %ilosci mozliwosci
    
    %% znajdywanie homozygoty
    
        haplotype = []; 
        freq = []; 
        for i = 1 : size(unphased,1)
            snp1 = unphased(i,1);
            snp2 = unphased(i,2);
            snp3 = unphased(i,3);
            if snp1 == 0 && snp2 == 0 && snp3 == 0
                idx = i;
                haplotype(1,:)= [3 2 3]; % {'G', 'C', 'G'};
                tf = ismember(unphased,[0 0 0],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            elseif snp1 == 2 && snp2 == 2 && snp3 == 2
                idx = i;
                haplotype(1,:)= [1 4 1]; % {'A', 'T', 'A'};
                tf = ismember(unphased,[2 2 2],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            elseif snp1 == 0 && snp2 == 0 && snp3 == 2
                idx = i;
                haplotype(1,:)= [3 2 1]; % {'G', 'C', 'A'};
                tf = ismember(unphased,[0 0 2],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            elseif snp1 == 0 && snp2 == 2 && snp3 == 0
                idx = i;
                haplotype(1,:)= [3 4 3]; % {'G', 'T', 'G'};
                tf = ismember(unphased,[0 2 0],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            elseif snp1 == 2 && snp2 == 0 && snp3 == 0
                idx = i;
                haplotype(1,:)= [1 2 3]; % {'A', 'C', 'G'};
                tf = ismember(unphased,[2 0 0],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            elseif snp1 == 0 && snp2 == 2 && snp3 == 2
                idx = i;
                haplotype(1,:)= [3 4 1]; % {'G', 'T', 'A'};
                tf = ismember(unphased,[0 2 2],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            elseif snp1 == 2 && snp2 == 0 && snp3 == 2
                idx = i;
                haplotype(1,:)= [1 2 1]; % {'A', 'c', 'A'};
                tf = ismember(unphased,[2 0 2],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            elseif snp1 == 2 && snp2 == 2 && snp3 == 0
                idx = i;
                haplotype(1,:)= [1 4 3]; % {'A', 'T', 'G'};
                tf = ismember(unphased,[2 2 0],'rows');
                freq(1) = patientsAmount(tf == 1) * 2;
                break;
            end
        end
        % pozybcie sie pierwszej homozygoty z danych
        unphased(i,:) = [];
        
        skip = true;
        %krok 3 szukanie nowego znanego itd.
        for i = 1:size(unphased,1)
           if (~skip)
           %sprawdzenia 
           end
           
           tmp = {};
           tmp = translate(unphased(6,:));
               
           if(size(tmp,1) == 1)
               homozygota
           else
               heterozygota
           end
           skip = false;
           
               
               
               
               
           
           
            
            
        end
        
        
    for i = 1:size(unphased,1)
        knownHaplo = translate(unphased(i,:), knownHaplo);
        knownHaplo{2,i} = patientsAmount(i);
    
    end
    
    
    %find homozygotes
    healthy0 = find(healthy(:,1) == 0 & healthy(:,2) == 0 & healthy(:,3) == 0);
    healthy2 = find(healthy(:,1) == 2 & healthy(:,2) == 2 & healthy(:,3) == 2);

    % label known haplotypes & count'em
    knownHaplotypes = healthy(healthy0(1),1:end);
    knownHaplotypes = [knownHaplotypes;healthy(healthy2(1),1:end)];
    [knownHaplo] = translate(knownHaplotypes(1,:), knownHaplo);
    knownHaplo{2,1} = size(healthy0,1);
    [knownHaplo] = translate(knownHaplotypes(2,:), knownHaplo);
    knownHaplo{2,2} = size(healthy2,1);

    for row = 1:size(healthy,1)
        tempRowHealthy = healthy(row,:);

        % tlumaczenie na zasady
        for i = 1:size(knownHaplotypes,1)
            [knownHaplo] = translate(knownHaplotypes(i,:), knownHaplo);
        end

        newHaploVariants = [];
        [newHaploVariants] = translate(tempRowHealthy, newHaploVariants);

        if size(unresolved,1) == 0
            unresolved = [unresolved; tempRowHealthy];
            tf = 0;
        else
            [tf] = ismember(unresolved,tempRowHealthy, 'rows');
        end

        found = [];
        for i = 1:size(knownHaplo,2)
            found = [found;strcmp(knownHaplo(1,i),newHaploVariants)];
        end
        bool = any(found);

        if (~any(bool) & tf)
            unresolved = [unresolved; tempRowHealthy]
        elseif (any(bool))
            idx = find(healthy(:,1) == tempRowHealthy(1,1) &...
                healthy(:,2) == tempRowHealthy(1,2) & ...
                healthy(:,3) == tempRowHealthy(1,3));
            amount = numel(idx);


            foundOne = newHaploVariants(bool == 1);

            for j = 1:size(foundOne,2)
                
                index = strcmp(knownHaplo(1,:),foundOne(1,j));

                %nie wiem jak to inaczej znalezc - rozwi¹zanie mega debilne juz
                %nawet nie glupie xDDD
                licznik = 1;
                for i = 1:size(index,2)
                    if(index(1,i) == 1) break;  end;
                    licznik = licznik + 1;
                end

                knownHaplo{2,licznik} = knownHaplo{2,licznik} + amount;

                newOnes = newHaploVariants(bool == 0);

                for i = 1:size(newOnes,2)
                    knownHaplo{1,end+1} = newOnes{:,i};
                    knownHaplo{2,end} = amount;
                end
            end
            %%wywalenie z unresolved todo
            idxUnresolved = find(unresolved(:,1) == tempRowHealthy(1,1) &...
                unresolved(:,2) == tempRowHealthy(1,2) & ...
                unresolved(:,3) == tempRowHealthy(1,3))
            unresolved(idxUnresolved,:) = [];
        else
            continue;
        end
    end
    out = knownHaplo;
end
