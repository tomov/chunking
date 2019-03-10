function plot_D(D)
    G = graph(D.G.E + D.G.hidden_E);
    h = plot(G, 'layout', 'force');

    for i = 1:length(D.G.hidden_edges)
        u = D.G.hidden_edges(i,1);
        v = D.G.hidden_edges(i,2);
        highlight(h, u, v, 'linestyle', '--');
    end
end
