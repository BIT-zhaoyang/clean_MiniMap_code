function D = connectGraph(dist)
% 'd' is the input genrated by function fast_knn.
% The generated graph may have disconnected components.
% This function checks if these disconnected components exit.
% If yes, connect them.

% The idea to connect these components is:
% 1. Find if there are any clusters by computing graphallshortestpaths
% 2. If there are Inf then it mean there are disconnected components.
% 3. Connect these components by a large distance.

% The function is partially copied from nystrom.m written by me.
%% ======================== Begin of function =============================
D = dist;
GDL = graphallshortestpaths(dist, 'directed', false);
if any(find(GDL == Inf))
       % If GDL has Inf elements, it means there are disconnected
       % components. 
%     disp('There are disconnected components');
    % Following three lines are used to find the number of components
    [tmp, firsts] = min(GDL == Inf);
    [comps, I, J] = unique(firsts);
    n_comps = length(comps);
    
    % Store different clusters in a cell structure
    Clusters = cell(n_comps,1);
    parfor i = 1:n_comps
        comps_index = find(firsts==comps(i));
        Clusters{i} = comps_index;                                     
    end  
    
    % Here is the connection implmentation. Since I only have 2 clusters
    % data, I can't guarantee it will successfully work on more than 2
    % clusters. But, it should be OK. You can try this.
    % The connection works like this. Say we have 3 clusters a, b, c. This
    % connection strategy connects a with b, and b with c. You can easily
    % modify it to make more connections between a and c. Just add another
    % for loop will be fine.
    for i = 1:n_comps - 1
        for j = i+1:n_comps
            EDM = dist(Clusters{i},Clusters{j}); % EDM: Euclidean distance 
                                                        % matrix between landmarks

            % Find the farthest landmarks in different clusters
            [farthest,ind] = max(EDM(:));
            [m,n] = ind2sub(size(EDM),ind);

            % Pay attention here. It may be a little tricky. In line 110, m,n
            % are index in EMD matrix. But, when we are going to connect
            % landmarks, we are operating in neighbourhood distance matrix. So,
            % we have to convert m,n to landmarks' real index.
            m = Clusters{i}(m);
            n = Clusters{j}(n);
            % Connect the pair landmarks with euclidean distance
            D(m,n) = max(max(dist));
        end
    end
    
    D = max(D,D');
    D = graphallshortestpaths(D);
    D = min(D, D');
else
    D = min(GDL, GDL');
end