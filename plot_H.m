function plot_H(H, D)
    G = graph(D.G.E);
    h = plot(G);

    cm = jet(max(H.c));
    for i = 1:D.G.N
        highlight(h, i, 'NodeColor', cm(H.c(i),:), 'MarkerSize', 10);
    end
end
