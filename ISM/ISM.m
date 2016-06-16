function Y = ISM(X, K, choice)
%% Compute neighbour
if ~exist('K', 'var')
    K = 5;
end

if ~exist('choice', 'var') | isempty(choice)
    choice = 'min';
end
[IDX, D] = fast_knn(X, K, choice);

%% Check disconnected components
D = connectGraph(D);

%% Sammon mapping
Y = mdscale(D, 2, 'criterion', 'sammon', 'start', 'random');
end