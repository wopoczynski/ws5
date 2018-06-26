<<<<<<< HEAD
function [chi2, D, chromosomCorelation] = X2(pq, P)
    D = zeros(2,3);
    chromosomCorelation = zeros(2,3);
    for i = 1:2
       for j = 1:3
          D(i,j) = P(j,1,i)*P(j,4,i) - P(j,2,i)*P(j,3,i);
          chromosomCorelation(i,j) = D(i,j)/ sqrt(pq(j,1,i)*pq(j,2,i)*pq(j,3,i)*pq(j,4,i));
        end
    end
    chi2 = chromosomCorelation .* 2;
end
=======
function [chi2, D, chromosomCorelation] = X2(pq, P)
    D = zeros(2,3);
    chromosomCorelation = zeros(2,3);
    for i = 1:2
       for j = 1:3
          D(i,j) = P(j,1,i)*P(j,4,i) - P(j,2,i)*P(j,3,i);
          chromosomCorelation(i,j) = D(i,j)/ sqrt(pq(j,1,i)*pq(j,2,i)*pq(j,3,i)*pq(j,4,i));
        end
    end
    chi2 = chromosomCorelation .* 2;
end
>>>>>>> ce0a84082a8cd9317b2111f7cfcd7d5378d33e4e
