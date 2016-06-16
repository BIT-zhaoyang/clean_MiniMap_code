function total = count_correct_neighbour(Y, label, K)
n = size(Y, 1);
err = zeros(n, 1);
[IDX, ~] = fast_knn(Y, K, 'max');
parfor i = 1:n
    neighbour = find(IDX(i, :));
    neighbour = label(neighbour);
    error = find(neighbour ~= label(i));
    err(i) = size(error, 1);
end
total = sum(err);
total = total / n / K * 100;
end