function [out] = clarkAlg(sample)

    knownHaplo = {};
    unresolved = [];

    tmp = sample(sample(:,1) == 1,:);
    healty = tmp(:,3:end);
    tmp = sample(sample(:,1) == 0,:);
    sic = tmp(:,3:end);

    %find homozygotes
    healty0 = find(healty(:,1) == 0 & healty(:,2) == 0 & healty(:,3) == 0);
    healty2 = find(healty(:,1) == 2 & healty(:,2) == 2 & healty(:,3) == 2);

    % label known haplotypes & count'em
    knownHaplotypes = healty(healty0(1),1:end);
    knownHaplotypes = [knownHaplotypes;healty(healty2(1),1:end)];
    [knownHaplo] = translate(knownHaplotypes(1,:), knownHaplo);
    knownHaplo{2,1} = size(healty0,1);
    [knownHaplo] = translate(knownHaplotypes(2,:), knownHaplo);
    knownHaplo{2,2} = size(healty2,1);

    for row = 1:size(healty,1)
        tempRowHealty = healty(row,:);

        % tlumaczenie na zasady
        for i = 1:size(knownHaplotypes,1)
            [knownHaplo] = translate(knownHaplotypes(i,:), knownHaplo);
        end

        newHaploVariants = [];
        [newHaploVariants] = translate(tempRowHealty, newHaploVariants);

        if size(unresolved,1) == 0
            unresolved = [unresolved; tempRowHealty];
            tf = 0;
        else
            [tf] = ismember(unresolved,tempRowHealty, 'rows');
        end

        found = [];
        for i = 1:size(knownHaplo,2)
            found = [found;strcmp(knownHaplo(1,i),newHaploVariants)];
        end
        bool = any(found);

        if (~any(bool) & tf)
            unresolved = [unresolved; tempRowHealty]
        elseif (any(bool))
            idx = find(healty(:,1) == tempRowHealty(1,1) &...
                healty(:,2) == tempRowHealty(1,2) & ...
                healty(:,3) == tempRowHealty(1,3));
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
            idxUnresolved = find(unresolved(:,1) == tempRowHealty(1,1) &...
                unresolved(:,2) == tempRowHealty(1,2) & ...
                unresolved(:,3) == tempRowHealty(1,3))
            unresolved(idxUnresolved,:) = [];
        else
            continue;
        end
    end
    out = knownHaplo;
end
