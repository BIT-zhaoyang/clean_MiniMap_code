function embed = nystrom(dist, P)
% function embed = nystrom(X, P,K)
%  This function takes Euclidean pairwise distance matrix 'dist' as input. 
%  Then, it construct the partial geodesic distance D.
%  The partition of the matrix takes form showing below.
%
%
%  Some notation for your convenience:
%  D: geodesic distance  
%  dist: pair-wise eudlidean distance   
%  dist_copy: as the literal mean
%
%  Nystr??m method description:
%      ___________        ___________
%      | A |  B  |        | E |  F  |
%      -----------        -----------
%  K = |   |     |    D = |   |     |
%      | B'|  C  |        | F'|  G  |
%      |___|_____|        |___|_____|
%  In the above graph, D is the complete geodesic distance. However, in
%  this function, D will be partial distace function, which is:
%       _____
%       | E |
%       -----
%   D = |   |
%       | F'|
%       |___|
%  The notation is derived from J.Platt(2005).
%  NOTICE that at first D is not squared. So I will reassign D to its
%  square later. I will still denote it as D.

%% 0. Set the number of landmarks
N = size(dist,1);   % Number of observations.
L = max(ceil(N * P),10);    % P represents the percentage of landmarks in the whole
                           % data.

D = zeros(L,N);     % D: geodesic distances
P = randperm(N);
idx = P(1:L);       % index of landmark points
IDX = 1:N;
IDX(idx) = [];      % index of the rest points

%% 1.construct neighbour graph
% find K nearest neighbour
dist_copy = dist;
%% 2. Connect different components if there are
% 1. First, compute the geodesic distance matrix between the L different
% landmark points.
% 2. If there are Inf, then it means there are differnet components.
% Compute how many different components there are. Then, connect the 
% farthest landmark points belong to different components.
% 3. After connecting different landmarks in different components. Compute
% the D*L submatrix. (This is implemented in section 3.)

% graphshortestpath(G,??S,??T) determines the single source-single 
% destination shortest path from node S to node T.
GDL = zeros(L,L); % GDL: geodesic distance between landmarks
parfor i = 1:L
    [d, path, pred] = graphshortestpath(dist,idx(i),idx);
    GDL(i,:) = d;
end
% Find if there are landmarks in different clusters
if any(find(GDL == Inf))
       % If GDL has Inf elements, it means there are disconnected
       % components. It may have the probability that all landmark points
       % sampled are in the same cluster so this method would fail. But
       % this is a very samll probability event so that we can ignore this.
       % If you are really unconfident about this, make the landmark points
       % more will be fine.
       
    disp('There are disconnected components');
    % Following three lines are used to find the number of components
    [tmp, firsts] = min(GDL == Inf);
    [comps, I, J] = unique(firsts);
    n_comps = length(comps);
    
    % Store landmarks in different clusters in a cell structure
    Clusters = cell(n_comps,1);
    parfor i = 1:n_comps
        comps_index = find(firsts==comps(i));
        Clusters{i} = idx(comps_index); % Here, we stroe the real index of
                                        % those landmarks in cell. By real 
                                        % index, I mean the index of landmarks
                                        % in given data. 
                                        % Pay attention to this!
    end
    
    
    % Here is the connection implmentation. Since I only have 2 clusters
    % data, I can't guarantee it will successfully work on more than 2
    % clusters. But, it should be OK. You can try this.
    % The connection works like this. Say we have 3 clusters a, b, c. This
    % connection strategy connects a with b, and b with c. You can easily
    % modify it to make more connections between a and c. Just add another
    % for loop will be fine.
    for i = 1:n_comps - 1
        EDM = dist_copy(Clusters{i},Clusters{i+1}); % EDM: Euclidean distance 
                                                    % matrix between landmarks
        
        % Find the farthest landmarks in different clusters
        [farthest,ind] = max(EDM(:));
        [m,n] = ind2sub(size(EDM),ind);
        
        % Pay attention here. It may be a little tricky. In line 110, m,n
        % are index in EMD matrix. But, when we are going to connect
        % landmarks, we are operating in neighbourhood distance matrix. So,
        % we have to convert m,n to landmarks' real index.
        m = Clusters{i}(m);
        n = Clusters{i+1}(n);
        % Connect the pair landmarks with euclidean distance
%         dist(m,n) = dist_copy(m,n);
        dist(m,n) = max(max(dist));
    end
    
    % Make the geodesic distance symmetric. It is the same as line 50.
    dist = max(dist,dist');
end
t = toc;
disp(['Computing partial geodesic distance uses ' num2str(t) ' seconds.'])
%% 3. Compute partial geodesic distance
% [dist, path, pred]??=??graphshortestpath(G,??S) determines the
% single-source shortest paths from node S to
% all other nodes in the graph represented by matrix G.
% Input G is an N-by-N sparse matrix that
% represents a graph. Nonzero entries in matrix G represent
% the weights of the edg51es. dist are the N distances
% from the source to every node (using Infs for nonreachable
% nodes and 0 for the source node). path contains
% the winning paths to every node. pred contains
% the predecessor nodes of the winning paths.
tic;
parfor i = 1:L
    [d, path, pred] = graphshortestpath(dist,idx(i));
    D(i,:) = d;
end
t = toc;
disp(['Computing partial geodesic distance uses ' num2str(t) ' seconds.'])
%% 4. Square the distance matrix
disp('Doing MDS...');
tic
D = D .^ 2; 

%% 5. Partition the matrix
E = D(1:L,idx);
F = D(1:L, IDX);

A = zeros(L,L);
B = zeros(L, N - L);

%% 6. Get A and B
A = double_centering(E);
A = 0.5 * (A + A'); % This step is to increase the stability of eigenvalue
                    % computing. Because, when I was testing the code. The
                    % eigenvalue of A may be complex eigenvalue. But a
                    % symmetric matrix shouldn't have complex eigenvalue. I
                    % suspect this is incured by doing double centering,
                    % MATLAB will lose some accuracy.
                    % So here, I force A to be symmetrical.
B = -0.5 * (F - repmat(sum(E,2)/L,1,N-L)); % Here, I used the method 
                                           % implemented by LMDS. Refer to
                                           % section 2.3 in J.Platt(2005).
                                           
%% 6*. How does original Isomap code handle the problem mentioned above:
% h = real(diag(val)); 
% [foo,sorth] = sort(h);  sorth = sorth(end:-1:1); 
% val = real(diag(val(sorth,sorth))); 
% vec = vec(:,sorth); 
% As we can see, the author only take the real part of the diagnol elments.
% So it seems he also encountered the same problem.

%% 7. Low dimensional embedding
K = 15;  % Number of dimensions to be preserved
[U,GAMA] = eigs(A, K,'LA');
% K = length(find(GAMA > 0.01)); % I am not sure whether it's again because of 
                            % the loss of accuracy. But MATLAB seems also
                            % perform not as well as I expected. So, I only
                            % retain eigenvalues which are bigger than 0.01
                            % rather than 0. 

A_embed = zeros(L,K);
V = zeros(L,K);
for i = 1:K
    A_embed(:,i) = sqrt(GAMA(i,i)) * U(:,i);
    V(:,i) = U(:,i) / sqrt(GAMA(i,i));
end

B_embed = B' * V;

embed = [A_embed; B_embed];
index = [idx IDX];
clear idx IDX
[useless, idx] = sort(index);
embed = embed(idx,:);
t = toc;
disp(['Doing MDS uses ' num2str(t) ' seconds.'])
end

%% 6**. Double centering function
function R = double_centering(E)
L = size(E,1);

e = ones(L,1);
I = eye(L);
H = I - e*e'/L;

R = - 0.5 * H * E * H;
end