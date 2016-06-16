clear, clc, close all;
addpath isomap tsne ISM MiniMap;

datanames = ReadLines('datalist.txt');
nd = length(datanames);

methods = { 'Isomap', ...
            'tSNE', ...
            'MiniMap'
          };
nm = length(methods);

set(gcf,'units','points','position',[0,0,800,800]);
left = 0.005;
bottom = 0.005;
width = 0.33;
height = 0.33;

K = 5; n = 5; choice = 'min';

for di = 1:nd
    dataname = datanames{di};
    fprintf('di = %d \t%10s\n', di, dataname);
    
    clear X label;
    [X, label] = GetData(dataname);
    if isempty(X)
        fprintf('No such data!\n');
        continue;
    end
    b = bottom + (di-1) * height;
    for mi = 1:nm
        method = methods{mi};
        fprintf('%s: ', method);
        switch method
            case 'Isomap'
                fprintf('\n');
                Y = Isomap(X, K, 'max');  % Isomap(X, K, choice)
            case 'I+SM'
                fprintf('\n');
                Y =  ISM(X, K, choice); % ISM(X, K, choice)
            case 'MiniMap'
                fprintf('\n');
                if(strcmp(dataname, 'USPS_partial.mat'))
                    K = 10;
                end
                Y = MiniMap(X, K, 'min', n);    % IsmDt(X, K, choice, nsteps, exagFactor, exagStep)
                K = 5;
            case 'tSNE'
                fprintf('\n');
                Y = fast_tsne(X, 2);
        end
        l = left + (mi - 1) * width;
        ax = axes('Position', [l b width height]);
        gscatter(Y(:, 1), Y(:, 2), label, [], [], 10, 'off');
        result = count_correct_neighbour(Y, label, 7);
        result = 100 - result;
        result = num2str(result,4)
        result = [result(1:4) '%'];
        text(0.80, 0.9, result, 'Units', 'normalized', 'FontSize', 14);
        set(ax, 'YTick', []);
        set(ax, 'XTick', []);
        xlabel('')
        ylabel('')     
    end
end