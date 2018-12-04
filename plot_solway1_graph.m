function h = plot_solway1_graph(H, D)

x = [-3 -1 -2 -3 -1 1 3 2 1 3];
y = [-3 -3 -2 -1 -1 1 1 2 3 3];
y = -y;
H.c = [1 1 1 1 1 2 2 2 2 2];


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
    else
        highlight(h, i, 'nodecolor', [247 148 30]/255);
    end

    for j = i+1:D.G.N
        if D.G.E(i,j)
            highlight(h, i, j, 'edgecolor', [0.1 0.1 0.1], 'linewidth', 2);
        end
    end
end

set(gca, 'xtick', []);
set(gca, 'ytick', []);
