function XcSquare = HW(counts)
rows = sum(counts(1,:,1));
p = (2 * counts(:,1,1) + counts(:,2,1) + counts(:,3,1)) ./ (rows *2);

q =  (2 * counts(:,4,1) + counts(:,2,1) + counts(:,3,1) ) ./ (rows *2);

p2 = p.^2;
pq2 = 2.*p.*q;
q2 = q.^2;

AA = counts(:,1,1);
Aa = counts(:,2,1) + counts(:,3,1);
aa = counts(:,4,1);

AAOcz = rows * p2;
AaOcz = rows * pq2;
aaOcz = rows * q2;

XcSquare = ((AA-AAOcz).^2 ./ AAOcz) + ((Aa-AaOcz).^2 ./ AaOcz ) + ((aa-aaOcz).^2 ./ aaOcz );
