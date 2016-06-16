function Y = IsmDt(X, K, choice, nsteps, exagFactor, exagStep)
% X: original data
%% Compute neighbour
% Normalize input data
X = X - min(X(:));
X = X / max(X(:));
X = bsxfun(@minus, X, mean(X, 1));

[n, m] = size(X);
if m > 50
    [~, project] = pca(X);
    X = project(:, 1:50);
end

if ~exist('K', 'var') | isempty(K)
    K = 5;
end

if ~exist('choice', 'var') | isempty(choice)
    choice = 'max';
end

% Selecting nearest neighbour is the most critical part in the algorithm
% here we used the idea in tSNE to select neighbours. This can be seen as
% selecting neighbours based on density.
[P, D] = x2p(X, K);
IDX = probNb(P, K);
%% Compute certain steps of geodesic distance
% Though the distance between one point and its neighbours should 
% be different, from results shown in tSNE, replace them by 1 should
% also work. If this is true, it should reveals that setting large
% distance to very large one is the critical point for tSNE to work.
if ~exist('nsteps', 'var') | isempty(nsteps)
    nsteps = 4;
end
GD = computeNeighbour(IDX, nsteps);
GD = max(GD, GD');
%% Set zero elements to be very large value
%  This corresponds to set small distances to constant small value
%  while set large distances to large constant value
n1 = sum(sum(GD == 1));
idx = (GD ~= 0);
% This part adjust the value for small distance
if exist('exagFactor', 'var') & ~isempty(exagFactor)
    GD(idx) = 1 / exagFactor;
else
    GD(idx) = log10(n)*log10(n) / n ;
end

n2 = sum(sum(GD == 0));
GD(GD == 0) = 1;
GD(logical(eye(size(GD)))) = 0;
GD = max(GD, GD');

%% Sammon mapping/MDS
opts.MaxIter = 10;
opts.Display = 'iter';
opts.Update = 'reduced';

% exageration for early training
w = GD .^ (-2); %.* n1 ./ n2 <- could be used when n1, n2 are very different
w(GD == 0) = 0; % set diagnol weights to be zero.
Y = smacof(GD, randn(size(X, 1), 2), w, opts);

% normal steps
w = GD .^ (-1); %.* n1 ./ n2;
w(GD == 0) = 0;
Y = smacof(GD, Y, w, opts);
end