function D = computeNeighbour(A, k)
% k represents the number of random walks.
if k == 1
    D = A;
    return;
end

D = A ^ k;
d = computeNeighbour(A, k - 1);

D = D + d;
D(find(D)) = 1;

D(logical(eye(size(D)))) = 0;
% D = max(D, D');
end