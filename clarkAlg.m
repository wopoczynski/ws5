function [haplotype, freq] = clarkAlg(sample)

% losowoœæ + zliczanie pacjentów
[unphased,~,ic] = unique(sample, 'rows');
patientsAmount = accumarray(ic, 1);
result = [unphased patientsAmount];
n = size(result,1);
idx = randperm(n);
result = result(idx,:);

unphased = result(:,1:end-1);
patientsAmount = result(:,end);

% krok 2 szukanie pierwszej homozygoty
haplotype = [];
freq = [];
for i = 1 : size(unphased,1)
    snp1 = unphased(i,1);
    snp2 = unphased(i,2);
    snp3 = unphased(i,3);
    if snp1 == 0 && snp2 == 0 && snp3 == 0
        haplotype(1,:)= [3 2 3];
        tf = ismember(unphased,[0 0 0],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    elseif snp1 == 2 && snp2 == 2 && snp3 == 2
        haplotype(1,:)= [1 4 1];
        tf = ismember(unphased,[2 2 2],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    elseif snp1 == 0 && snp2 == 0 && snp3 == 2
        haplotype(1,:)= [3 2 1];
        tf = ismember(unphased,[0 0 2],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    elseif snp1 == 0 && snp2 == 2 && snp3 == 0
        haplotype(1,:)= [3 4 3];
        tf = ismember(unphased,[0 2 0],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    elseif snp1 == 2 && snp2 == 0 && snp3 == 0
        haplotype(1,:)= [1 2 3];
        tf = ismember(unphased,[2 0 0],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    elseif snp1 == 0 && snp2 == 2 && snp3 == 2
        haplotype(1,:)= [3 4 1];
        tf = ismember(unphased,[0 2 2],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    elseif snp1 == 2 && snp2 == 2 && snp3 == 0
        haplotype(1,:)= [1 4 3];
        tf = ismember(unphased,[2 2 0],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    elseif snp1 == 2 && snp2 == 0 && snp3 == 2
        haplotype(1,:)= [1 2 1];
        tf = ismember(unphased,[2 0 2],'rows');
        freq(1) = patientsAmount(tf == 1) * 2;
        break;
    end
end

unphased(i,:) = [];
% krok 3 - szukanie nowych, zapamietanie  brakujacych itd.
keep = [];
nonFirst = 0;
for i = 1:size(unphased,1)
    if nonFirst == 1
        lengthKeep = size(keep, 2);
        for q = 1:lengthKeep
            tmpKeep = cell2mat(keep(q));
            lengthTmp = size(tmpKeep,1);
            if lengthTmp > 0
                lenH = size(haplotype,1);
                for k = 1:lenH
                    for j = 1:lengthTmp
                        tfTmp(j,k) = isequal(tmpKeep(j,:), haplotype(k,:));
                    end
                end
                if sum(sum(tfTmp)) == 1
                    if lengthTmp == 1
                        freq(sum(tfTmp,1) == 1) = freq(sum(tfTmp,1) == 1) + patientsAmount(i) * 2;
                        keep{1,q} = keep{1,q} * 5;
                    else % HETEROZYGOTA
                        lengthTmp = size(tmp,1);
                        lengthHaplotype = size(haplotype,1);
                        lengthTF = size(tfOld,2);
                        if lengthTF == 1 % homozygota
                            if sum(sum(tfOld,2) == 0 ) == 1
                                haplotype(lengthHaplotype+1, :) = tmp(sum(tfOld,2) == 0,:);
                                freq(lengthHaplotype+1) =  patientsAmount(i);
                            else
                                idx = sum(tfOld,1) == 1;
                                freq(idx) = freq(idx) + patientsAmount(i);
                            end
                        elseif lengthTF > 1 % heterozygota
                            appended = 0;
                            for q = 1:2:lengthTmp
                                if (sum(sum(tfOld(q:q+1,:),2)) == 1) && (appended == 0)
                                    idx = sum(tfOld(q:q+1,:),2) == 0;
                                    t = tmp(q:q+1,:);
                                    for k = 1:size(haplotype,1)
                                        tf3(k) = isequal(t(idx,:), haplotype(k,:));
                                    end
                                    if sum(tf3) == 0
                                        haplotype(size(haplotype,1)+1, :) = t(idx,:);
                                        freq(size(haplotype,1)+1) =  patientsAmount(i);
                                    else
                                    end
                                    idx = sum(tfOld(q:q+1,:),1);
                                    idx = find(idx == 1);
                                    freq(idx) = freq(idx) + patientsAmount(i);
                                    appended = 1;
                                elseif sum(sum(tfOld(q:q+1,:),2)) == 2
                                    idx = sum(tfOld(q:q+1,:),1);
                                    idx = find(idx == 1);
                                    freq(idx) = freq(idx) + patientsAmount(i);
                                    
                                end
                            end
                        end
                        keep{1,q} = keep{1,q} * 5;
                    end
                end
                tmpKeep = [];
                tfOld = [];
            else
            end
        end
    end
    tmp = [];
    tmp = translate(unphased(i,:));
    if size(tmp,1) == 1 % HOMOZYGOTA
        for k = 1:size(haplotype,1)
            tf(k) = isequal(tmp, haplotype(k,:));
        end
        if sum(tf) == 0
            keep{1,i} = tmp;
            freq(tf == 1) = freq(tf == 1) + patientsAmount(i) * 2;
        end
    else % HETEROZYGOTA
        [haplotype, freq, keep, tf] = hetero(tmp, haplotype, freq, keep, i, patientsAmount);
        tfOld = tf;
    end
    nonFirst = 1;
end