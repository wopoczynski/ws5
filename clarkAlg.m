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
        tempRowHealty = healty(row,:); %%przeklejone usunac

        % tlumaczenie na zasady
        for i = 1:size(knownHaplotypes,1)
            [knownHaplo] = translate(knownHaplotypes(i,:), knownHaplo);
        end
        
        newHaploVariants = [];
        [newHaploVariants] = translate(tempRowHealty, newHaploVariants);
        
        if size(unresolved,1) == 0
            unresolved = [unresolved; tempRowHealty]
            tf = 0;
        else 
            [tf] = ismember(unresolved,tempRowHealty, 'rows'); 
        end
        
        if ((~any(strcmp(knownHaplo,newHaplo)) ||...
           ~any(strcmp(knownHaplo,newHaploA))) &...
            tf)
            unresolved = [unresolved; tempRowHealty]
        else
            continue;
        end
  
        tempRowHealty = healty(row,:); %%przeklejone usunac
        idx = find(healty(:,1) == tempRowHealty(1,1) &...
                healty(:,2) == tempRowHealty(1,2) & ...
                healty(:,3) == tempRowHealty(1,3));

        knownHaplotypesAmount = [size(healty0,1);size(healty2,1)];

    end
    
   
    
    
    
    
    out = 'todo';
end
