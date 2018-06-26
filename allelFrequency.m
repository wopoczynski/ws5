<<<<<<< HEAD
function [counts, frequency, pq, recombination, P] = allelFrequency(sample)
    counts = zeros(3,4,2);
    frequency = counts;
    pq = counts;
    P = counts;
    recombination = counts;
    for j = 1:2
        X = sample(sample(:,1) == (j-1) ,:);
        for i=1:3
            counts(i,1,j) = sum(X(:,i*2)==1 & X(:,(i*2+1))==1); % AB
            counts(i,2,j) = sum(X(:,i*2)==1 & X(:,(i*2+1))==2); % Ab
            counts(i,3,j) = sum(X(:,i*2)==2 & X(:,(i*2+1))==1); % aB
            counts(i,4,j) = sum(X(:,i*2)==2 & X(:,(i*2+1))==2); % ab
        end
        frequency(:,:,j) = counts(:,:,j) ./ 200;
        for i=1:3
            pq(i,1,j) = frequency(i,1,j) + frequency(i,2,j); % p1
            pq(i,2,j) = frequency(i,3,j) + frequency(i,4,j); % p2
            pq(i,3,j) = frequency(i,1,j) + frequency(i,3,j); % q1
            pq(i,4,j) = frequency(i,2,j) + frequency(i,4,j); % q2
            recombination(j,i) = frequency(i,2,j) + frequency(i,3,j);
            P(i,[1,4],j) = (1-recombination(j,i))/2; % p11, p22
            P(i,[2,3],j) = recombination(j,i)/2; % p12, p21
        end
    end
=======
function [counts, frequency, pq, recombination, P] = allelFrequency(sample)
    counts = zeros(3,4,2);
    frequency = counts;
    pq = counts;
    P = counts;
    recombination = counts;
    for j = 1:2
        X = sample(sample(:,1) == (j-1) ,:);
        for i=1:3
            counts(i,1,j) = sum(X(:,i*2)==1 & X(:,(i*2+1))==1); % AB
            counts(i,2,j) = sum(X(:,i*2)==1 & X(:,(i*2+1))==2); % Ab
            counts(i,3,j) = sum(X(:,i*2)==2 & X(:,(i*2+1))==1); % aB
            counts(i,4,j) = sum(X(:,i*2)==2 & X(:,(i*2+1))==2); % ab
        end
        frequency(:,:,j) = counts(:,:,j) ./ 200;
        for i=1:3
            pq(i,1,j) = frequency(i,1,j) + frequency(i,2,j); % p1
            pq(i,2,j) = frequency(i,3,j) + frequency(i,4,j); % p2
            pq(i,3,j) = frequency(i,1,j) + frequency(i,3,j); % q1
            pq(i,4,j) = frequency(i,2,j) + frequency(i,4,j); % q2
            recombination(j,i) = frequency(i,2,j) + frequency(i,3,j);
            P(i,[1,4],j) = (1-recombination(j,i))/2; % p11, p22
            P(i,[2,3],j) = recombination(j,i)/2; % p12, p21
        end
    end
>>>>>>> ce0a84082a8cd9317b2111f7cfcd7d5378d33e4e
end