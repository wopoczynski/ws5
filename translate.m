function tmp = translate(unphased)
    n = sum(unphased == 1);
    options = 2.^n;

    %% SNP 1 GG/AG/AA
    if unphased(1) == 0
        tmp(1:options,1) = 'G';
    elseif unphased(1) == 2
        tmp(1:options,1) = 'A';
    elseif unphased(1) == 1
        idxNParz = 1:2:options;
        idxParz = 2:2:options+1;
        tmp(idxNParz, 1) = 'A';
        tmp(idxParz, 1) = 'G';
    end
    %% SNP 2 CC/CT/TT
    if unphased(2) == 0
        tmp(1:options,2) = 'C';
    elseif unphased(2) == 2
        tmp(1:options,2) = 'T';
    elseif unphased(2) == 1
        if n == 1
            tmp(1, 2) = 'C';
            tmp(2, 2) = 'T';
        elseif n == 2
            tmp([1 4], 2) = 'C';
            tmp([2 3], 2) = 'T';
        else
            tmp([1 3 6 8], 2) = 'C';
            tmp([2 4 5 7], 2) = 'T';
        end
    end
    %% SNP3 GG/AG/AA
    if unphased(3) == 0
        tmp(1:options, 3) = 'G';
    elseif unphased(3) == 2
        tmp(1:options, 3) = 'A';
    elseif unphased(3) == 1
        if n == 1
            tmp(1, 3) = 'A';
            tmp(2, 3) = 'G';
        elseif n == 2
            tmp([1 4], 3) = 'A';
            tmp([2 3], 3) = 'G';
        else
            tmp([1 4 5 8], 3) = 'A';
            tmp([2 3 6 7], 3) = 'G';
        end
    end
end