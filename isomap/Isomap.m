function Y = Isomap(X, K, choice)
%% Compute neighbour
n = size(X, 1);
if nargin < 2 | isempty(K)
    K = 15;
end

if ~exist('choice', 'var') | isempty(choice)
    choice = 'max';
end
[IDX, D] = fast_knn(X, K, choice);

%% Check disconnected components
D = connectGraph(D);

%% Compute geodesic distance
% D = graphallshortestpaths(D, 'Directed', 'false');
% D = min(D, D');
%% MDS
Y = cmdscale(D);
end