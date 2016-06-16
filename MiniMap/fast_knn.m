function [ind, dist] = fast_knn(X, k, choice, verbose)
% Calculate approximated kNN in a fast way
% A Matlab wrapper that calls van der Maaten's ball tree implementation
% 
%   knn = fast_approx_knn(X, k)
%   knn = fast_approx_knn(X, k, verbose)
%
% Input:
%   X : N x d
%   k : the number of neighbors
%   verbose : to display the progress or not, default=false
% Output:
%   knn : N x N
%
% Copyright (c) 2014, Zhirong Yang (Aalto University)
% All rights reserved.

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

n = size(X,1);

X = X';

J = reshape(repmat(1:n, k, 1), n*k, 1);
[I, D] = fastknn_mex(X, k, verbose);
ind = sparse(I, J, ones(n*k,1), n, n, n*k)';
dist = sparse(I, J, D, n, n, n*k)';

% I will modify ind & dist here. 
% 26.01.2015    Yang Zhao
if strcmp(choice, 'max')
    ind = max(ind, ind');
    dist = max(dist, dist');
end
if strcmp(choice, 'min')
    ind = min(ind, ind');
    dist = min(dist, dist');
end
end