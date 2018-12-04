clear all;
load('solway4_alpha=2_10000.mat');

figure('pos', [100 100 1000 600] * 3/4);
fontsize = 13;
axisfontsize = 10;
lettersize = 20;

% A: graph
%
subplot(2,3,1);

[h, xs, ys] = plot_solway4_graph(H, D);
labelnode(h, 10, ' start');
labelnode(h, 20, ' goal');
title('Experimental Design', 'fontsize', fontsize);

hold on;
E = zeros(size(D.G.E));
E(10,5) = 1;
E(5,6) = 1;
E(6,8) = 1;
E(8,9) = 1;
E(9,19) = 1;
E(19,20) = 1;
G = digraph(E);
h = plot(G);
set(h, 'XData', xs);
set(h, 'YData', ys);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    highlight(h, i, 'marker', 'none');
    for j = 1:D.G.N
        if E(i,j)
            highlight(h, i, j, 'edgecolor', [244 137 143]/255, 'linewidth', 2);
        end
    end
end
h.ArrowSize = 10;
h.EdgeAlpha = 1;
%h.ArrowPosition = 0.9



E = zeros(size(D.G.E));
E(10,12) = 1;
E(12,16) = 1;
E(16,18) = 1;
E(18,23) = 1;
E(23,22) = 1;
E(22,20) = 1;
G = digraph(E);
h = plot(G);
set(h, 'XData', xs);
set(h, 'YData', ys);
labelnode(h, 1:D.G.N, '');
for i = 1:D.G.N
    highlight(h, i, 'marker', 'none');
    for j = 1:D.G.N
        if E(i,j)
            highlight(h, i, j, 'edgecolor', [143 207 165]/255, 'linewidth', 2);
        end
    end
end
h.ArrowSize = 10;
h.EdgeAlpha = 1;
%h.ArrowPosition = 0.9

hold off;




% B: Data
%

subplot(2,3,2);

n = size(tasks,1) * N; 
c1 = 0.72 * n; 
c2 = n - c1;
p = 2 * binocdf(min(c1,c2), n, 0.5);
y = binoinv([0.025 0.975], n, 0.5) / n;

m = c1/n;
bar(m);
hold on;
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [y(1) y(1) y(2) y(2)], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');

set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
xticklabels({'P(fewer boundaries)'});
hold off;

title('Data', 'fontsize', fontsize);


% C: Model

subplot(2,3,3);

c1 = 0;
c2 = 0;
n = 0;

synth = [];
for t = 1:size(tasks,1)
    c1 = c1 + sum(move(:,t) == nexts(t,1)); % count 1
    c2 = c2 + sum(move(:,t) == nexts(t,2)); % count 2
    synth = [synth ones(1,c1) zeros(1,c2)];
end
n = c1 + c2;
p = 2 * binocdf(min(c1,c2), n, 0.5);
y = binoinv([0.025 0.975], n, 0.5) / n;

m = c1/n;
se = sem(synth);
bar(m);
hold on;
errorbar(m, se);
line([0 2], [0.5 0.5], 'linestyle', '--', 'color', [0.6 0.6 0.6]);
h = fill([0 2 2 0], [y(1) y(1) y(2) y(2)], [0.4 0.4 0.4]);
set(h, 'facealpha', 0.5, 'edgecolor', 'none');

set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [0 1]);
set(gca, 'ytick', [0 0.5 1]);
xticklabels({'P(fewer boundaries)'});
hold off;

title('Model', 'fontsize', fontsize);


% D: Hierarchies

x = xs;
y = ys;
for s = 1:12
    subplot(4,6, 12 + s);
    h = plot_H(map_H{s}, D);

    set(h, 'XData', x);
    set(h, 'YData', y);
    labelnode(h, 1:D.G.N, '');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    h.MarkerSize = 6;

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
print('solway4.pdf', '-dpdf');

% stats
%
fprintf('two-sided binomial test n = %d, p = %.4f\n', n, p);
