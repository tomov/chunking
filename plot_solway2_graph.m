function h = plot_solway2_graph(H, D)

x = [-3 -2 -1 -3 -2 -1 -3 -2 -1 0 1 2 3 1 2 3 1 2 3];
y = [-1 -1 -1 0 0 0 1 1 1 0 -1 -1 -1 0 0 0 1 1 1];
H.c = [1 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 3 3 3];


h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    highlight(h, i, 'nodecolor', [0.1 0.1 0.1], 'markersize', 12);
end
hold on;

h = plot_H(H, D);
set(h, 'XData', x);
set(h, 'YData', y);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    if H.c(i) == 1
        highlight(h, i, 'nodecolor', [16 118 188]/255);
    elseif H.c(i) == 3
        highlight(h, i, 'nodecolor', [247 148 30]/255);
    else
        highlight(h, i, 'nodecolor', [113 84 159]/255);
    end

    for j = i+1:D.G.N
        if D.G.E(i,j)
            highlight(h, i, j, 'edgecolor', [0.1 0.1 0.1], 'linewidth', 2);
        end
    end
end

set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gca, 'ylim', [-4 2]);
set(gca, 'xlim', [-4 4]);
