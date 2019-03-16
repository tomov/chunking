
clear all;


figure('pos', [10 500 600 1100] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;



plots_n = 10;
plots_m = 6;
plot_idx = 0;

%
% schapiro
%

load schapiro_N=30_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat  % <-- sample_c

r1 = 6;
r2 = 3;
cx = cos([0:2]*2*pi/3 + pi/2) * r1;
cy = sin([0:2]*2*pi/3 + pi/2) * r1;
x(1:5) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(1);
y(1:5) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(1);
x(6:10) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(2);
y(6:10) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(2);
x(11:15) = cos([0:4]*2*pi/5 + pi/2) * r2 + cx(3);
y(11:15) = sin([0:4]*2*pi/5 + pi/2) * r2 + cy(3);


for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);
    h = plot_H(chosen_H{s}, D);

    set(h, 'XData', x);
    set(h, 'YData', y);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    h.MarkerSize = 4;

    if s == 3
        title(['                            ', 'Schapiro et al. (2013)'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end

%
% solway 1
%

load solway1_N=40_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat  % <-- sample_c

x = [-3 -1 -2 -3 -1 1 3 2 1 3];
y = [-3 -3 -2 -1 -1 1 1 2 3 3];
y = -y;
for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);
    h = plot_H(chosen_H{s}, D);

    set(h, 'XData', x);
    set(h, 'YData', y);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    h.MarkerSize = 4;

    if s == 3
        title(['                            ', 'Solway et al. (2014), Experiment 1'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end

%
% solway 2
%

load solway2_N=10_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat % sample_c

x = [-3 -2 -1 -3 -2 -1 -3 -2 -1 0 1 2 3 1 2 3 1 2 3];
y = [-1 -1 -1 0 0 0 1 1 1 0 -1 -1 -1 0 0 0 1 1 1];
for s = 1:10
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);
    h = plot_H(chosen_H{s}, D);

    set(h, 'XData', x);
    set(h, 'YData', y);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'ylim', [-3 3]);
    h.MarkerSize = 4;

    if s == 3
        title(['                            ', 'Solway et al. (2014), Experiment 2'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end

plot_idx = plot_idx + 2; % only 10 subjects here...


%
% solway 4
%

load solway4_N=35_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat  % <-- sample_c

for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);
    [h, xs, ys] = plot_solway4_graph(H, D);
    hold off;
    h = plot_H(chosen_H{s}, D);

    set(h, 'XData', xs);
    set(h, 'YData', ys);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    h.MarkerSize = 4;

    if s == 3
        title(['                            ', 'Solway et al. (2014), Experiment 4'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end


%
% lynn
%

load lynn_N=78_alpha=1.0000_nsamples=10000_eps=0.6000_last.mat  % <-- sample_c

for s = 1:12
    plot_idx = plot_idx + 1;
    subplot(plots_n, plots_m, plot_idx);
    h = plot_H(chosen_H{s}, D);

    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    h.MarkerSize = 4;

    if s == 3
        title(['                            ', 'Lynn et al. (2019)'], 'fontsize', fontsize);
    end
    xlabel(sprintf('subject %d', s));
end




% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'portrait');
print('figures/hierarchies_Sim.pdf', '-dpdf');
