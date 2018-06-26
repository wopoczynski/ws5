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
            
            idx = find(knownHaplo(1,:)
            %% todo 
            % w bool znajduje siê true/false elementu odnalezionego w
            % znanych haplotypach, trzeba wyjac znaleziony dodac now¹ iloœæ
            % wyst¹pieñ do neigo, dopisaæ pozosta³e i do pozosta³ych iloœæ
            % wyst¹pieñ, a po wszystkich usun¹æ z unresolved i next
            % iteracja

            %%wywalenie z unresolved todo
            idxUnresolved = find(unresolved(:,:) == tempRowHealty(:,:))
            unresolved(idx,:) = [];
        else
            continue;
        end
        
       

    end
    
   
    
    
    
    
    out = 'todo';
end
