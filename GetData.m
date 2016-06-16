function [X, label] = GetData(dataname)
filename = ['./data/', dataname];
if isempty(dir(filename))
    X = [];
    label = [];
    return;
end

load(filename);

n = size(label, 1);
idx = randperm(n);

label = label(idx,:);
X = X(idx,:);
end