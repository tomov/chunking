clear all;
%load('solway2.mat');
%load('solway2_alpha=5_nsamples=1000.mat'); % <-- preprint
load solway2_N=10_alpha=5.0000_nsamples=1000.mat  % <-- sample repros!!!!
%load('solway2_alpha=2_nsamples=100.mat');

figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

% A: graph
%
subplot(2,3,1);

h = plot_solway2_graph(H, D);
d = 0.1;
text(h.XData(9), h.YData(9)  , 's', 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(h.XData(16) , h.YData(16) + d, 'g', 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%labelnode(h, 9, ' start');
%labelnode(h, 16, ' goal');
title('Experimental Design', 'fontsize', fontsize);


% B: Data
%

sizes = (100 - 74) * ones(1,19) / 18;
sizes(10) = 74;
sizes = sizes * 0.15;

ax = subplot(2,3,2);
h = plot_solway2_graph(H, D);
for i = 1:D.G.N
    highlight(h, i, 'NodeColor', [0.6 0.6 0.6], 'MarkerSize', 12 + sizes(i));
end

hold on;
plot_solway2_graph(H, D);
hold off;

title('Data', 'fontsize', fontsize);


axes('Position', [ax.Position(1) ax.Position(2) 0.07 0.12]);
box on;

for j = 1:null_iters
    null_p(j) = mean(null{j}(corr(:)) == 10);
end
null_p = sort(null_p);
lcb = null_p(length(null_p) * 0.025);
ucb = null_p(length(null_p) * 0.975);

m = 0.74;

bar(m);
hold on;
line([0 2], [mean(null_p) mean(null_p)], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [lcb lcb ucb ucb], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
xticklabels({'P(choose bottleneck)'});
hold off;



% C: Model

ax = subplot(2,3,3);
h = plot_solway2_graph(H, D);
for i = 1:D.G.N
    f = sum(loc(corr(:)) == i);
    f = f / length(corr(:)) * sum(sizes);
    highlight(h, i, 'NodeColor', [0.6 0.6 0.6], 'MarkerSize', 10 + f);
end

hold on;
plot_solway2_graph(H, D);
hold off;

title('Model', 'fontsize', fontsize);


axes('Position', [ax.Position(1) ax.Position(2) 0.07 0.12]);
box on;

for j = 1:null_iters
    null_p(j) = mean(null{j}(corr(:)) == 10);
end
null_p = sort(null_p);
lcb = null_p(length(null_p) * 0.025);
ucb = null_p(length(null_p) * 0.975);

m = mean(p);
se = sem(p);

bar(m);
hold on;
errorbar(m, se, 'color', 'black');
line([0 2], [mean(null_p) mean(null_p)], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
%h = fill([0 2 2 0], [lcb lcb ucb ucb], [0.4 0.4 0.4]);
%set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
xticklabels({'P(choose bottleneck)'});
hold off;




% D: Hierarchies

x = [-3 -2 -1 -3 -2 -1 -3 -2 -1 0 1 2 3 1 2 3 1 2 3];
y = [-1 -1 -1 0 0 0 1 1 1 0 -1 -1 -1 0 0 0 1 1 1];
for s = 1:10
    subplot(4,5, 10 + s);
    h = plot_H(map_H{s}, D);

    set(h, 'XData', x);
    set(h, 'YData', y);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'ylim', [-3 3]);

    if s == 3
        title('Example hierarchies', 'fontsize', fontsize);
    end
end


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1);
text(0.10, 0.96, 'A', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.36, 0.96, 'B', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.64, 0.96, 'C', 'FontSize', lettersize, 'FontWeight', 'bold');
text(0.10, 0.52, 'D', 'FontSize', lettersize, 'FontWeight', 'bold');


% save figure
h = gcf;
%set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
%print('figures/solway2.pdf', '-dpdf');


% stats
%
pos = find(m <= null_p);
if isempty(pos) 
    pos = 0;
end
pos = pos(end);
fprintf('%f; MC test (%d samples from null), p = %e\n', m, null_iters, pos / length(null_p));
