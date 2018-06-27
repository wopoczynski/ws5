function [haplotype, freq, keep, tf] = hetero(tmp, haplotype, freq, keep, i, patientsAmount)

for a = 1:size(haplotype,1)
    for j = 1:size(tmp,1)
        tf(j,a) = isequal(tmp(j,:), haplotype(a,:));
    end
end
if sum(sum(tf)) == 0
    keep{1,i} = tmp;
elseif sum(sum(tf)) > 0
    lengthTmp = size(tmp,1);
    lengthHaplotype = size(haplotype,1);
    lengthTF = size(tf,2);
    if lengthTF == 1 % homozygota
        if sum(sum(tf,2) == 0 ) == 1
            haplotype(lengthHaplotype+1, :) = tmp(sum(tf,2) == 0,:);
            freq(lengthHaplotype+1) =  patientsAmount(i);
        else
            idx = sum(tf,1) == 1;
            freq(idx) = freq(idx) + patientsAmount(i);
        end
    elseif lengthTF > 1 % heterozygota
        appended = 0;
        for b = 1:2:lengthTmp
            if (sum(sum(tf(b:b+1,:),2)) == 1) && (appended == 0)
                idx = sum(tf(b:b+1,:),2) == 0;
                t = tmp(b:b+1,:);
                for a = 1:size(haplotype,1)
                    tf3(a) = isequal(t(idx,:), haplotype(a,:));
                end
                if sum(tf3) == 0
                    haplotype(size(haplotype,1)+1, :) = t(idx,:);
                    freq(size(haplotype,1)+1) =  patientsAmount(i);
                else
                end
                idx = sum(tf(b:b+1,:),1);
                idx = find(idx == 1);
                freq(idx) = freq(idx) + patientsAmount(i);
                appended = 1;
            elseif sum(sum(tf(b:b+1,:),2)) == 2
                idx = sum(tf(b:b+1,:),1);
                idx = find(idx == 1);
                freq(idx) = freq(idx) + patientsAmount(i);
                
            end
        end
    end
end
end