
figure('pos', [10 500 900 690]);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;



plots_n = 7;
plots_m = 12;
plot_idx = 0;

%
% subway 10 (from fig_subway10.m)
%

load('model_Exp_1_thru_4_samples=10000_alpha=1.0000_eps=0.6000_last.mat');
i = 3;
j = 1;

for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);

    H = pl(i).H{j}(s);
    D = pl(i).D{j}(s);

    [~, xs, ys] = plot_subway10_graph(H, D); % just to get xs & ys
    hold off;
    h = plot_H(H, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-1 3]);
    h.MarkerSize = 6;

    if s == 6
        title(['                            ', 'Experiment 1'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end

%
% subway 9 (fig_subway9.m)
%

i = 4;
j = 1:3;

j = [1 1 1 1 1 1 3 3 3 3 3 3];
for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);

    H = pl(i).H{j(s)}(s);
    D = pl(i).D{j(s)}(s);
    [~, xs, ys] = plot_subway9_graph(H, D); % just to get xs & ys
    hold off;
    h = plot_H(H, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-1 3]);
    h.MarkerSize = 6;

    if s == 6
        title(['                            ', 'Experiment 2'], 'fontsize', fontsize);
    end
    if s == 1
        ylabel('bad');
    end
    if s == 7
        ylabel('good');
    end
    xlabel(sprintf('subject %d', s));
end

%
% unlearn (fig_unlearn.m)
%

load('model_exp_v2_3_circ_alpha=1.0000_nsamples=10000_div_eps=0.6000_last.mat');

idx = [1 2 3 4 5 6 1 2 3 4 5 6];
subj = [1 1 1 1 1 1 2 2 2 2 2 2];
for i = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);

    s = subj(i);

    H = chosen_H{s, idx(i)};
    D = D_full(s);
    [h, xs, ys] = plot_unlearn_graph(H, D);
    hold off;
    h = plot_H(H, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-1 3]);
    h.MarkerSize = 6;

    if i == 6
        title(['                            ', 'Experiment 3'], 'fontsize', fontsize);
    end
    if i == 1
        ylabel('subject 1');
    end
    if i == 7
        ylabel('subject 2');
        xlabel('probe #1');
    end
    if i >= 8
        xlabel(sprintf('probe #%d', i - 6));
    else
        xlabel(sprintf('probe #%d', i));
    end
end


%
% subway 10 map (fig_subway10_map.m)
%

load('model_Exp_1_thru_4_samples=10000_alpha=1.0000_eps=0.6000_last.mat');

i = 1;
j = 1;
for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);

    H = pl(i).H{j}(s);
    D = pl(i).D{j}(s);
    [h, xs, ys] = plot_subway10_graph(H, D);
    hold off;
    h = plot_H(H, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-1 3]);
    h.MarkerSize = 6;

    if s == 6
        title(['                            ', 'Experiment 4'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end


%
% subway 9 map (fig_subway9_map.m)
%

i = 2;
j = [1 1 1 2 2 2 3 3 3 2 2 2];
for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);

    H = pl(i).H{j(s)}(s);
    D = pl(i).D{j(s)}(s);
    [h, xs, ys] = plot_subway9_graph(H, D);
    hold off;
    h = plot_H(H, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-1 3]);
    h.MarkerSize = 6;

    if s == 6
        title(['                            ', 'Experiment 5'], 'fontsize', fontsize);
    end
    if s == 1
        ylabel('bad');
    end
    if s == 4
        ylabel('control 1');
    end
    if s == 7
        ylabel('control 2');
    end
    if s == 10
        ylabel('good');
    end
    xlabel(sprintf('subject %d', s));
end


%
% mines (fig_mines.m)
%


load('mines_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat'); % sample_c

for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);

    [h, xs, ys] = plot_mines_graph(H, D);
    hold off;
    h = plot_H(chosen_H{s}, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'ylim', [-3 3]);
    h.MarkerSize = 6;

    if s == 6
        title(['                            ', 'Experiment 6'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end


%
% mines 10 map (fig_mines10_map.m
%

load('mines10_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat');

i = 5;
j = 1;

for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);

    [h, xs, ys] = plot_subway10_graph(H, D);
    hold off;
    h = plot_H(chosen_H{s}, D);
    set(h, 'XData', xs);
    set(h, 'YData', ys);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlim', [-1 3]);
    h.MarkerSize = 6;

    if s == 6
        title(['                            ', 'Experiment 7'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end



% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
print('figures/hierarchies_Exp.pdf', '-dpdf');
