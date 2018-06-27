function [haplotype, freq] = clarkAlg(sample)

    %% krok 1 zliczanie wystêpuj¹cych i losowanie kolejnosci
    [unphased,~,ic] = unique(sample, 'rows');
    patientsAmount = accumarray(ic, 1);
    result = [unphased patientsAmount];
    n = size(result,1);
    idx = randperm(n);
    result = result(idx,:);

    unphased = result(:,1:end-1); % mozliwosci
    patientsAmount = result(:,end); %ilosci mozliwosci

    % znajdywanie homozygoty
    haplotype = {};
    freq = [];
    for i = 1 : size(unphased,1)
        snp1 = unphased(i,1);
        snp2 = unphased(i,2);
        snp3 = unphased(i,3);
        if snp1 == 0 && snp2 == 0 && snp3 == 0
            idx = i;
            haplotype(1,:) = {'GCG'};
            tf = ismember(unphased,[0 0 0],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        elseif snp1 == 2 && snp2 == 2 && snp3 == 2
            idx = i;
            haplotype(1,:) = {'ATA'};
            tf = ismember(unphased,[2 2 2],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        elseif snp1 == 0 && snp2 == 0 && snp3 == 2
            idx = i;
            haplotype(1,:) = {'GCA'};
            tf = ismember(unphased,[0 0 2],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        elseif snp1 == 0 && snp2 == 2 && snp3 == 0
            idx = i;
            haplotype(1,:) = {'GTG'};
            tf = ismember(unphased,[0 2 0],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        elseif snp1 == 2 && snp2 == 0 && snp3 == 0
            idx = i;
            haplotype(1,:) = {'ACG'};
            tf = ismember(unphased,[2 0 0],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        elseif snp1 == 0 && snp2 == 2 && snp3 == 2
            idx = i;
            haplotype(1,:) = {'GTA'};
            tf = ismember(unphased,[0 2 2],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        elseif snp1 == 2 && snp2 == 0 && snp3 == 2
            idx = i;
            haplotype(1,:) = {'AcA'};
            tf = ismember(unphased,[2 0 2],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        elseif snp1 == 2 && snp2 == 2 && snp3 == 0
            idx = i;
            haplotype(1,:) = {'ATG'};
            tf = ismember(unphased,[2 2 0],'rows');
            freq(1) = patientsAmount(tf == 1) * 2;
            break;
        end
    end

    % pozybcie sie pierwszej homozygoty z danych
    unphased(i,:) = [];
    
    tf = [];
    freq = [];
    keep = [];
    skip = true;
    %% krok 3 szukanie nowego znanego itd.
    for i = 1:size(unphased,1)
        if (~skip)
            dimKeep = size(keep, 2);
            for j = 1:dimKeep
                %% sprawdzenie jakiœ haplotyp znajduje siê w macierzy hap
                x = cell2mat(keep(j));
                dimX = size(x,1);
                if dimX > 0
                    haplotypeSize = size(haplotype,1);
                    for k = 1:haplotypeSize
                        for l = 1:dimX
                            tfTemp(l,k) = isequal(x(l,:), haplotype(k,:));
                        end
                    end
                    if any(any(tfTemp)) == 1 
                        if dimX == 1 % HOMOZYGOTA 
                            freq(sum(tfTemp,1) == 1) = freq(sum(tfTemp,1) == 1) + patientsAmount(i) * 2;
                            keep{1,j} = keep{1,j} * 5;
                        else 
                            [haplotype, freq] = sprawdzanie_dopisywanie (tf, temp, haplotype, freq, i, patientsAmount);
                            keep{1,j} = keep{1,j} * 5;
                        end
                    end
                    x = [];
                    tf = [];
                end    
            end
        end
        
        tmp = {};
        tmp = translate(unphased(5,:));
        if(size(tmp,1) == 1)            %homozygota
            tf = ismember(tmp, char(haplotype), 'rows');
            if (any(tf) == 0)
                keep{1,i} = tmp;
            else
                 freq(tf == 1) = freq(tf == 1) + patientsAmount(i) * 2;
            end
        else %heterozygota
            for j = 1:size(haplotype,1)
                tf(:,j) = ismember(tmp, char(haplotype), 'rows');
            end
            if (any(any(tf)) == 0)
                keep{1,i} = tmp;
            else
                if (size(tf,2) == 1) %homozygota nowa
                    if any(any(tf) == 0) == 1
                        haplotypeTmp = char(haplotype);
                        haplotypeTmp(size(haplotypeTmp,1)+1, :) = tmp(sum(tf,2) == 0,:);
                        freq(size(haplotypeTmp,1)+1) =  patientsAmount(i);
                    else
                        idx = sum(tf,1) == 1;
                        freq(idx) = freq(idx) + patientsAmount(i);
                    end
                else % heterozygota nowa
                    stop = false; % warunkuje zatrzymanie dopisywania haplotypów po dopisaniu jednego
                    for j = 1:2:size(tmp,1)
                        tf(j:j+1,:)
                        if (sum(sum(tf(j:j+1,:),2)) == 1) && (stop == false)
                            idx = sum(tf(j:j+1,:),2) == 0;
                            t = tmp(j:j+1,:);
                            % weryfikacja czy na pewno nowy
                            for k = 1:size(haplotype,1)
                                tfTmp(k) = isequal(t(idx,:), haplotype(k,:));
                            end
                            if sum(tfTmp) == 0
                                haplotypeTmp = char(haplotype);
                                haplotypeTmp(size(haplotypeTmp,1)+1, :) = t(idx,:);
                                freq(size(haplotype,1)+1) =  patientsAmount(i);
                            else
                            end
                            % sumowanie zliczeñ
                            idx = sum(tf(j:j+1,:),1);
                            idx = find(idx == 1);
                            freq(idx) = freq(idx) + patientsAmount(i);
                            stop = true;
                        elseif any(any(tf(j:j+1,:),2)) == 2
                            % sumowanie zliczeñ dla ka¿dego z 2
                            idx = sum(tf(j:j+1,:),1);
                            idx = find(idx == 1);
                            freq(idx) = freq(idx) + patientsAmount(i);
                        end
                    end
                end
            end
        skip = false;
        end    
    end
end
