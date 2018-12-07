clear all;
%load('solway1.mat');
load('solway1_alpha=2_nsamples=100.mat');

figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

% A: graph
%
subplot(2,3,1);

plot_solway1_graph(H, D);
title('Experimental Design', 'fontsize', fontsize);


% B: Data
%

sizes = [3 3 3 2 10 12 3 1 2 1] * 2;

ax = subplot(2,3,2);
h = plot_solway1_graph(H, D);
for i = 1:D.G.N
    highlight(h, i, 'NodeColor', [0.6 0.6 0.6], 'MarkerSize', 12 + sizes(i));
end

hold on;
plot_solway1_graph(H, D);
hold off;

title('Data', 'fontsize', fontsize);


axes('Position', [ax.Position(1) ax.Position(2) 0.07 0.12]);
box on;

c1 = floor(length(loc(:)) * (4.4 / 5.4));
c2 = length(loc(:)) - c1;
n = c1 + c2;
p = 1 - binocdf(c1, n, 2/10);
ci = binoinv([0.025 0.975], n, 2/10) / n;

m = c1/n;
bar(m);
hold on;
line([0 2], [2/10 2/10], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [ci(1) ci(1) ci(2) ci(2)], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
xticklabels({'P(choose bottlenecks)'});

hold off;


% C: Model

ax = subplot(2,3,3);
h = plot_solway1_graph(H, D);
for i = 1:D.G.N
    f = sum(loc(:) == i);
    f = f / length(loc(:)) * sum(sizes);
    highlight(h, i, 'NodeColor', [0.6 0.6 0.6], 'MarkerSize', 10 + f);
end

hold on;
plot_solway1_graph(H, D);
hold off;

title('Model', 'fontsize', fontsize);

axes('Position', [ax.Position(1) ax.Position(2) 0.07 0.12]);
box on;

c1 = sum(loc(:) == 5 | loc(:) == 6); % count 1
c2 = sum(loc(:) ~= 5 & loc(:) ~= 6); % count 2
n = c1 + c2;
p = 1 - binocdf(c1, n, 2/10);
ci = binoinv([0.025 0.975], n, 2/10) / n;

m = c1/n;
se = std(loc(:) == 5 | loc(:) == 6) / sqrt(n);
bar(m);
hold on;
errorbar(m, se, 'color', 'black');
line([0 2], [2/10 2/10], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [ci(1) ci(1) ci(2) ci(2)], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');
set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
xticklabels({'P(choose bottlenecks)'});

hold off;



% D: Hierarchies

x = [-3 -1 -2 -3 -1 1 3 2 1 3];
y = [-3 -3 -2 -1 -1 1 1 2 3 3];
y = -y;
for s = 1:12
    subplot(4,6, 12 + s);
    h = plot_H(map_H{s}, D);

    set(h, 'XData', x);
    set(h, 'YData', y);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);

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
print('solway1.pdf', '-dpdf');

% stats
%
fprintf('right-tailed binomial test %f, n = %d, p = %e\n', c1/n, n, p);