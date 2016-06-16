function mat2txt(name)
% X: matrix     name: name of the dataset
filename = [name '.txt'];
dataname = ['./data/' name '.mat'];

% read data
load(dataname);

% write size of matrix
[m, n] = size(X);
fileID = fopen(filename, 'w');
fprintf(fileID, '%d %d\n', m, n);
fclose(fileID);

% write data from X to txt
dlmwrite(filename, X, '-append', 'delimiter', ' ', 'precision', '%.4f');

end