function h = plot_Hs(H, D)

    n = ceil(length(H) / 10);
    m = round(length(H) / n);

    figure;

    for l = 1:length(H)

        subplot(n, m, l);

        G = graph(D.G.E + D.G.hidden_E);
        h = plot(G, 'layout', 'force');

        for i = 1:size(D.G.hidden_edges,1)
            u = D.G.hidden_edges(i,1);
            v = D.G.hidden_edges(i,2);
            highlight(h, u, v, 'linestyle', '--');
        end

        cm = jet(max(H(l).c));
        for i = 1:D.G.N
            highlight(h, i, 'NodeColor', cm(H(l).c(i),:), 'MarkerSize', 5);
        end
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
end
