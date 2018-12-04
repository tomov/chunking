function [h] = plot_lynn_graph(H, D)

H.c = [1 2 2 3 3 4 4 4 4 4 4 3 3 2 2];


h = plot_H(H, D);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    highlight(h, i, 'nodecolor', [0.1 0.1 0.1], 'markersize', 12);
end
hold on;

h = plot_H(H, D);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    if H.c(i) == 1
        highlight(h, i, 'nodecolor', [255 255 255]/255);
    elseif H.c(i) == 2
        highlight(h, i, 'nodecolor', [29 117 0]/255);
    elseif H.c(i) == 3
        highlight(h, i, 'nodecolor', [0 118 186]/255);
    else
        highlight(h, i, 'nodecolor', [181 23 0]/255);
    end

    for j = i+1:D.G.N
        if D.G.E(i,j)
            highlight(h, i, j, 'edgecolor', [0.1 0.1 0.1], 'linewidth', 2);
        end
    end
end

hold off;

set(gca, 'xtick', []);
set(gca, 'ytick', []);
