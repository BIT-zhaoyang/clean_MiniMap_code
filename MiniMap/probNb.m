function A = probNb(P, k)
%% This function abandons some points if they are not close neighbours to a given point
% D is distance matrix corresponding to the adjacency matrix.
% A is the new adjacency matrix after probability processing
% We call the selected neighbour after this process the core neighbour.

ny = size(P, 1);
A = zeros(ny, ny);

parfor i = 1:ny
    [ele, I] = sort(P(i, :), 'descend');
    for j = 3:ny
        if ele(j) <= (0.2 * ele(1))
            row = zeros(1, ny);
            row(I(1:j - 1)) = 1;
            A(i, :) = row;
            break;
        end
    end
end

end